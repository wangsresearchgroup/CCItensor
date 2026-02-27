required_pkgs <- c("reshape2", "ggplot2", "igraph", "rTensor", "ggrepel", "umap")
maybe_install(required_pkgs)
quiet_library(required_pkgs)

# ----------------------------- Args ----------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_tensor_pipeline.R <input_rds> <output_dir> <component>", call. = FALSE)
}

input_rds <- args[[1]]
out_dir   <- args[[2]]
component <- args[[3]]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "pipeline.log")
zz <- file(log_file, open = "wt")
sink(zz, split = TRUE)
on.exit({
  try(sink(), silent = TRUE)
  try(close(zz), silent = TRUE)
}, add = TRUE)

cat("Starting pipeline\n")
cat("Input     :", input_rds, "\n")
cat("Output dir:", out_dir, "\n")
cat("Component :", component, "\n\n")
set.seed(100)

# -------------------------- Helper functions --------------------------
stop_if_missing_cols <- function(df, cols) {
  missing <- setdiff(cols, colnames(df))
  if (length(missing) > 0) {
    stop("Missing columns in table: ", paste(missing, collapse = ", "))
  }
}

# Frequency of (source,target,ligand,receptor) across samples
find_frequency <- function(table_list) {
  keys <- unlist(lapply(table_list, function(tab) {
    if (is.null(tab) || nrow(tab) == 0) return(character(0))
    stop_if_missing_cols(tab, c("source", "target", "ligand", "receptor"))
    paste(tab$source, tab$target, tab$ligand, tab$receptor, sep = "_")
  }), use.names = FALSE)
  
  table(keys)
}

filter_top_n_keys <- function(table_list, freq_table, N = 1000) {
  if (length(freq_table) == 0) return(table_list)
  keep_keys <- names(sort(freq_table, decreasing = TRUE))[seq_len(min(N, length(freq_table)))]
  
  lapply(table_list, function(tab) {
    if (is.null(tab) || nrow(tab) == 0) return(tab)
    stop_if_missing_cols(tab, c("source", "target", "ligand", "receptor"))
    keys <- paste(tab$source, tab$target, tab$ligand, tab$receptor, sep = "_")
    tab[keys %in% keep_keys, , drop = FALSE]
  })
}

five_tensor <- function(input_list) {
  # Convert list of per-sample CCI tables to a 5D array:
  # source x target x ligand x receptor x sample, values = prob
  if (length(input_list) == 0) stop("input_list is empty")
  sample_names <- names(input_list)
  if (is.null(sample_names) || any(sample_names == "")) {
    stop("input_list must be a named list of samples.")
  }
  
  # Combine
  combined <- do.call(rbind, lapply(input_list, as.data.frame))
  stop_if_missing_cols(combined, c("source", "target", "ligand", "receptor", "prob"))
  
  sources   <- unique(combined$source)
  targets   <- unique(combined$target)
  ligands   <- unique(combined$ligand)
  receptors <- unique(combined$receptor)
  
  arr <- array(
    0,
    dim = c(length(sources), length(targets), length(ligands), length(receptors), length(input_list)),
    dimnames = list(sources, targets, ligands, receptors, sample_names)
  )
  
  # Fill per sample (vectorized using matrix indexing)
  for (i in seq_along(input_list)) {
    tab <- as.data.frame(input_list[[i]])
    if (nrow(tab) == 0) next
    stop_if_missing_cols(tab, c("source", "target", "ligand", "receptor", "prob"))
    
    idx <- cbind(
      match(tab$source,   sources),
      match(tab$target,   targets),
      match(tab$ligand,   ligands),
      match(tab$receptor, receptors),
      i
    )
    arr[idx] <- tab$prob
  }
  
  list(tensor = arr, sources = sources, targets = targets, ligands = ligands, receptors = receptors, samples = sample_names)
}

construct_inner_matrix <- function(cp_res, tensor_meta) {
  # Inner = t( (SxL) %*% t(TxR) )
  U <- cp_res$U
  
  rownames(U[[1]]) <- tensor_meta$sources
  rownames(U[[2]]) <- tensor_meta$targets
  rownames(U[[3]]) <- tensor_meta$ligands
  rownames(U[[4]]) <- tensor_meta$receptors
  rownames(U[[5]]) <- tensor_meta$samples
  
  outer_by_component <- function(A, B) {
    lapply(seq_len(ncol(A)), function(k) outer(A[, k], B[, k]))
  }
  
  flatten_named <- function(mat) {
    v <- as.vector(mat)
    dn <- dimnames(mat)
    if (!is.null(dn) && length(dn) == 2) {
      names(v) <- as.vector(outer(dn[[1]], dn[[2]], paste, sep = "_"))
    }
    v
  }
  
  array_list_to_table <- function(arr_list) {
    do.call(cbind, lapply(arr_list, flatten_named))
  }
  
  s_l <- array_list_to_table(outer_by_component(U[[1]], U[[3]]))
  t_r <- array_list_to_table(outer_by_component(U[[2]], U[[4]]))
  
  t(s_l %*% t(t_r))
}

generate_top_signatures <- function(mat, frac = 0.0002) {
  stopifnot(frac > 0, frac <= 1)
  df <- data.frame(
    ligand = rep(colnames(mat), each = nrow(mat)),
    receptor = rep(rownames(mat), times = ncol(mat)),
    value = as.vector(mat)
  )
  df <- df[order(df$value, decreasing = TRUE), , drop = FALSE]
  head(df, ceiling(nrow(df) * frac))
}

create_adjacency_matrix <- function(mat, threshold = 0.1) {
  # threshold is fraction of max(mat); outputs 0/1 matrix
  stopifnot(threshold >= 0, threshold <= 1)
  abs_thresh <- max(mat) * threshold
  (mat >= abs_thresh) * 1
}

make_vertex_colors <- function(vertex_names) {
  cell_types <- sub("_.*$", "", vertex_names)
  uniq <- unique(cell_types)
  
  # Simple palette you can edit; recycled if needed
  pal <- c(
    "#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a",
    "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6",
    "#b15928", "#ffff99", "#8dd3c7", "#bebada", "#fb8072"
  )
  
  cols <- setNames(rep(pal, length.out = length(uniq)), uniq)
  cols[cell_types]
}

create_graph <- function(inner_matrix, threshold_frac = 0.1) {
  adj <- create_adjacency_matrix(inner_matrix, threshold_frac)
  g <- igraph::graph_from_biadjacency_matrix(adj, mode = "in", directed = TRUE, weighted = NULL)
  
  # Remove isolates
  g <- igraph::delete_vertices(g, igraph::degree(g) == 0)
  if (igraph::vcount(g) == 0) return(g)
  
  vnames <- igraph::V(g)$name
  lr_names <- sub("^[^_]*_", "", vnames)
  vcols <- make_vertex_colors(vnames)
  
  # Use full inner_matrix weights for edges
  el <- igraph::as_data_frame(g, what = "edges")
  w <- mapply(function(from, to) inner_matrix[to, from], el$from, el$to)
  igraph::E(g)$weight <- w
  
  # Styling
  igraph::V(g)$label <- lr_names
  igraph::V(g)$label.cex <- 1.0
  igraph::V(g)$label.color <- vcols
  igraph::V(g)$color <- "lightyellow"
  
  # Scale node sizes by centrality computed on weighted full graph
  full_g <- igraph::graph_from_biadjacency_matrix(inner_matrix, mode = "in", directed = TRUE, weighted = TRUE)
  cen <- igraph::eigen_centrality(full_g)$vector
  cen <- cen[igraph::V(full_g)$name %in% igraph::V(g)$name]
  igraph::V(g)$size <- pmax(2, as.numeric(cen) * 40)
  
  g
}

run_umap_mode5 <- function(cp_res, n_neighbors = 7, k = 3) {
  set.seed(100)
  factors <- cp_res$U[[5]]
  um <- umap::umap(factors, n_neighbors = n_neighbors)
  emb <- um$layout
  
  set.seed(100)
  clusters <- kmeans(emb, centers = k)$cluster
  
  df <- data.frame(UMAP1 = emb[, 1], UMAP2 = emb[, 2], Cluster = factor(clusters))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(UMAP1, UMAP2, color = Cluster)) +
    ggplot2::geom_point(size = 3, alpha = 0.7) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right") +
    ggplot2::stat_ellipse(ggplot2::aes(group = Cluster), linetype = "dashed", color = "black")
  
  list(plot = p, clusters = clusters, embedding = emb)
}

# ---------------------------- Load data -------------------------------
all_components <- readRDS(input_rds)
if (!component %in% names(all_components)) {
  stop("component not found in input_rds. Available: ", paste(names(all_components), collapse = ", "))
}

table_list <- all_components[[component]]
sample_names <- names(table_list)

if (length(table_list) <= 1) stop("Need >1 sample tables in component list.")
cat("Samples:", length(table_list), "\n")

# ------------------------ Frequency filtering -------------------------
cat("Computing key frequencies...\n")
freq <- find_frequency(table_list)

N <- 1000
cat("Keeping top N keys:", N, "\n\n")
filtered_table_list <- filter_top_n_keys(table_list, freq, N = N)

# ------------------------ Tensor + decomposition -----------------------
cat("Constructing tensor...\n")
tensor_meta <- five_tensor(filtered_table_list)
tensor_5 <- rTensor::as.tensor(tensor_meta$tensor)

cat("Running CP decomposition...\n")
cp_res <- rTensor::cp(tensor_5, num_components = 20, max_iter = 50)
saveRDS(cp_res, file.path(out_dir, "cp_result.rds"))
rm(tensor_5); gc()

# ------------------------ Inner matrix + heatmap -----------------------
cat("Computing inner matrix...\n")
inner_matrix <- construct_inner_matrix(cp_res, tensor_meta)
saveRDS(inner_matrix, file.path(out_dir, "inner_matrix.rds"))

cat("Saving heatmap...\n")
p_heat <- ggplot2::ggplot(reshape2::melt(inner_matrix), ggplot2::aes(Var2, Var1, fill = value)) +
  ggplot2::geom_tile() +
  ggplot2::theme_minimal()
ggplot2::ggsave(file.path(out_dir, "inner_heatmap.png"), p_heat, width = 10, height = 8)

# ------------------------ Top signatures ------------------------------
cat("Saving top interactions...\n")
top_sigs <- generate_top_signatures(inner_matrix, frac = 0.0002)
write.csv(top_sigs, file.path(out_dir, "top_pos_interactions.csv"), row.names = FALSE)

# ------------------------ Network graph -------------------------------
cat("Saving network graph...\n")
g <- create_graph(inner_matrix, threshold_frac = 0.1)
png(file.path(out_dir, "network_graph.png"), width = 3100, height = 3100, res = 300)
plot(g)
dev.off()

# ------------------------ UMAP ---------------------------------------
cat("Saving UMAP...\n")
um <- run_umap_mode5(cp_res, n_neighbors = 7, k = 3)
png(file.path(out_dir, "umap.png"), width = 1500, height = 1500, res = 300)
print(um$plot)
dev.off()

cat("\nPipeline completed successfully\n")