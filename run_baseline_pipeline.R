required_pkgs <- c(
  "reshape2",
  "ggplot2",
  "igraph",
  "ggrepel",
  "umap"
)

installed <- rownames(installed.packages())
to_install <- setdiff(required_pkgs, installed)
if (length(to_install) > 0) {
  install.packages(to_install, dependencies = TRUE, quiet = TRUE)
}

suppressPackageStartupMessages({
  library(reshape2)
  library(ggplot2)
  library(igraph)
  library(ggrepel)
  library(umap)
})

# ---------- Args ----------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_baseline_pipeline.R <input_rds> <output_dir> <component>")
}

input_rds <- args[1]
out_dir   <- args[2]
component <- args[3]
meta_rds <- args[4]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
sink(file.path(out_dir, "pipeline.log"), split = TRUE)

cat("Baseline CCI pipeline\n")
cat("Input:", input_rds, "\n")
cat("Component:", component, "\n")

set.seed(100)

# ============================================================
# Load data
# ============================================================
meta <- readRDS(meta_rds)
all_components <- readRDS(input_rds)
stopifnot(component %in% names(all_components))

table_list <- all_components[[component]]
table_list <- table_list[!vapply(table_list, is.null, logical(1))]
stopifnot(length(table_list) > 1)

sample_names <- names(table_list)
print("sample names")
print(sample_names)
# ============================================================
# Helper: table → LR matrix
# rows = target_cell_receptor
# cols = source_cell_ligand
# ============================================================

make_lr_matrix <- function(tab, value_col = "prob") {
  stopifnot(all(c("source","target","ligand","receptor", value_col) %in% colnames(tab)))
  
  row_key <- paste(tab$target, tab$receptor, sep = "_")
  col_key <- paste(tab$source, tab$ligand, sep = "_")
  
  rows <- unique(row_key)
  cols <- unique(col_key)
  
  mat <- matrix(
    0,
    nrow = length(rows),
    ncol = length(cols),
    dimnames = list(rows, cols)
  )
  
  for (i in seq_len(nrow(tab))) {
    mat[row_key[i], col_key[i]] <- mat[row_key[i], col_key[i]] + tab[[value_col]][i]
  }
  mat
}

# ============================================================
# Build per-sample LR matrices
# ============================================================

lr_mats <- lapply(table_list, make_lr_matrix)

all_rows <- unique(unlist(lapply(lr_mats, rownames)))
all_cols <- unique(unlist(lapply(lr_mats, colnames)))

align_matrix <- function(mat) {
  out <- matrix(0, nrow = length(all_rows), ncol = length(all_cols),
                dimnames = list(all_rows, all_cols))
  out[rownames(mat), colnames(mat)] <- mat
  out
}

lr_mats <- lapply(lr_mats, align_matrix)

# ============================================================
# Baseline inner matrix = mean across samples
# ============================================================

cat("Computing mean LR interaction matrix\n")
inner_matrix <- Reduce(`+`, lr_mats) / length(lr_mats)
saveRDS(inner_matrix, file.path(out_dir, "inner_matrix_no_tensor.rds"))

# ============================================================
# Heatmap
# ============================================================

p_heat <- ggplot(melt(inner_matrix), aes(Var2, Var1, fill = value)) +
  geom_tile() +
  theme_minimal() +
  labs(title = paste("Mean LR matrix (baseline)", component))

ggsave(
  file.path(out_dir, "inner_heatmap_no_tensor.png"),
  p_heat, width = 10, height = 8
)

# ============================================================
# Network graph
# ============================================================

create_adjacency_matrix <- function(mat, frac = 0.1) {
  thresh <- max(mat) * frac
  (mat >= thresh) * 1
}

# Cell type = everything before first "_"
vertex_cell_types <- unique(
  gsub("_.*$", "", c(colnames(inner_matrix), rownames(inner_matrix)))
)

vertex_colors <- c(
  "yellowgreen","#67153d","darkgreen","#FF4500","#8B0000",
  "#4682B4","#9400D3","#6d4e3e","darkorange","pink","black",
  "#1f78b4","#33a02c","#e31a1c","#ff7f00","#6a3d9a",
  "#b15928","#a6cee3","#b2df8a","#fb9a99","#fdbf6f",
  "#cab2d6","#ffff99","#8dd3c7","#bebada","#fb8072",
  "#80b1d3","#fdb462","#b3de69","#fccde5"
)
vertex_colors <- rep(vertex_colors, length.out = length(vertex_cell_types))
names(vertex_colors) <- vertex_cell_types

adj <- create_adjacency_matrix(inner_matrix, 0.1)
g <- graph_from_biadjacency_matrix(adj, mode = "in", directed = TRUE)

V(g)$cell_type <- gsub("_.*$", "", V(g)$name)
V(g)$label <- gsub("^[^_]*_", "", V(g)$name)
V(g)$label.cex <- 1.2
V(g)$label.color <- vertex_colors[V(g)$cell_type]
V(g)$color <- "lightyellow"

centrality <- eigen_centrality(
  graph_from_biadjacency_matrix(inner_matrix, mode = "in", weighted = TRUE)
)$vector
V(g)$size <- ifelse(centrality > 0, centrality * 40, 1)

g <- delete_vertices(g, degree(g) == 0)

png(file.path(out_dir, "network_graph_no_tensor.png"),
    width = 3100, height = 3100, res = 300)
plot(g)
dev.off()

# ============================================================
# UMAP on vectorized LR matrices
# ============================================================
# ------------------------------------------------------------
# Normalize sample IDs to match metadata donor_id
# ------------------------------------------------------------
sample_ids <- gsub("_cellchat$", "", sample_names)

cat("Running baseline UMAP\n")

X <- do.call(rbind, lapply(lr_mats, as.vector))
rownames(X) <- sample_ids

print(sample_ids)
print(meta$donor_id)

set.seed(100)
umap_res <- umap(X, n_neighbors = 7)

umap_df <- data.frame(
  donor_id = sample_ids,
  UMAP1 = umap_res$layout[,1],
  UMAP2 = umap_res$layout[,2],
  stringsAsFactors = FALSE
)

# ---- Join metadata (FIXED)
umap_df <- umap_df |>
  dplyr::left_join(meta, by = "donor_id")

# ---- Sanity check (check a REAL metadata column, not donor_id)
if (any(is.na(umap_df$sex))) {
  bad <- umap_df$donor_id[is.na(umap_df$sex)]
  stop(
    "Metadata missing for some samples in baseline UMAP:\n",
    paste(head(bad, 10), collapse = ", ")
  )
}

# ---- Clustering
set.seed(100)
umap_df$cluster <- factor(kmeans(umap_res$layout, centers = 3)$cluster)
cluster_tbl <- data.frame(
  donor_id = umap_df$donor_id,
  cluster  = as.integer(umap_df$cluster),
  stringsAsFactors = FALSE
)
write.csv(
  cluster_tbl,
  file.path(out_dir, "cluster_membership.csv"),
  row.names = FALSE
)



# ---- Cluster UMAP
p_umap <- ggplot(umap_df, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_classic() +
  stat_ellipse(aes(group = cluster), linetype = "dashed") +
  labs(title = paste("UMAP on vectorized LR matrices (baseline)", component))

ggsave(
  file.path(out_dir, "umap_no_tensor.png"),
  p_umap, width = 6, height = 6
)

# ============================================================
# Metadata UMAPs (baseline, no tensor)
# ============================================================

plot_discrete <- function(df, col, title) {
  
  p <- ggplot(df, aes(UMAP1, UMAP2, color = .data[[col]])) +
    geom_point(size = 2.2, alpha = 0.75) +
    theme_classic() +
    ggtitle(title)
  
  # ---- Ellipses ONLY for valid groups
  df_ell <- df |>
    dplyr::mutate(
      ellipse_group = dplyr::case_when(
        col == "uicc_stage" & .data[[col]] %in% c("I","II","1","2") ~ "Stage 1/2",
        col == "uicc_stage" & .data[[col]] %in% c("III","IV","3","4") ~ "Stage 3/4",
        col != "uicc_stage" ~ as.character(.data[[col]]),
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::filter(
      !is.na(ellipse_group),
      !(ellipse_group %in% c("unknown","Unknown","NA"))
    )
  
  if (nrow(df_ell) > 0 && length(unique(df_ell$ellipse_group)) >= 2) {
    p <- p +
      stat_ellipse(
        data = df_ell,
        aes(x = UMAP1, y = UMAP2, fill = ellipse_group),
        geom = "polygon",
        alpha = 0.25,
        color = NA,
        inherit.aes = FALSE
      )
  }
  
  p
}

plot_continuous <- function(df, col, title) {
  ggplot(df, aes(UMAP1, UMAP2, color = .data[[col]])) +
    geom_point(size = 2.2, alpha = 0.75) +
    theme_classic() +
    ggtitle(title)
}

meta_discrete   <- c("sex", "uicc_stage", "ever_smoker")
meta_continuous <- c("age")

for (col in meta_discrete) {
  if (col %in% colnames(umap_df)) {
    ggsave(
      file.path(out_dir, paste0("umap_", col, "_no_tensor.png")),
      plot_discrete(
        umap_df,
        col,
        paste("Baseline UMAP:", component, "-", col)
      ),
      width = 6, height = 6
    )
  }
}

for (col in meta_continuous) {
  if (col %in% colnames(umap_df) && is.numeric(umap_df[[col]])) {
    ggsave(
      file.path(out_dir, paste0("umap_", col, "_no_tensor.png")),
      plot_continuous(
        umap_df,
        col,
        paste("Baseline UMAP:", component, "-", col)
      ),
      width = 6, height = 6
    )
  }
}


cat("Baseline pipeline complete\n")
sink()
