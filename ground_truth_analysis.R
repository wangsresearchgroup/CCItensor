SIM_DIR <- "simulations"
OUT_DIR <- "simulations/ground_truth"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# Load lig_params per sample
# -------------------------

sample_dirs <- list.dirs(
  SIM_DIR,
  recursive = FALSE,
  full.names = TRUE
)

sample_dirs <- sample_dirs[grepl("^sample_", basename(sample_dirs))]

if (length(sample_dirs) == 0) {
  stop("No sample_* directories found in ", SIM_DIR)
}

lig_params_list <- lapply(sample_dirs, function(sdir) {
  f <- file.path(sdir, "lig_params.csv")
  
  if (!file.exists(f)) {
    stop("Missing lig_params.csv in ", sdir)
  }
  
  df <- read.csv(f, stringsAsFactors = FALSE)
  
  required_cols <- c("target", "regulator", "effect")
  if (!all(required_cols %in% colnames(df))) {
    stop("lig_params.csv in ", sdir,
         " must contain columns: ", paste(required_cols, collapse = ", "))
  }
  
  df |>
    select(target, regulator, effect)
})

names(lig_params_list) <- basename(sample_dirs)

# Enforce stable ordering
lig_params_list <- lig_params_list[order(names(lig_params_list))]

# Save list
saveRDS(
  lig_params_list,
  file.path(OUT_DIR, "lig_params_list.rds")
)

cat("Saved lig_params_list.rds with",
    length(lig_params_list), "samples\n")

# -------------------------
# Assemble ground-truth summary
# -------------------------

GT_summary <- bind_rows(
  lapply(names(lig_params_list), function(s) {
    lig_params_list[[s]] |>
      mutate(sample = s)
  })
) |>
  group_by(target, regulator) |>
  summarise(
    freq        = n_distinct(sample),
    mean_effect = mean(effect),
    sum_effect  = sum(effect),
    .groups = "drop"
  ) |>
  arrange(desc(sum_effect))

# Save outputs
saveRDS(
  GT_summary,
  file.path(OUT_DIR, "ground_truth_summary.rds")
)

write.csv(
  GT_summary,
  file.path(OUT_DIR, "ground_truth_summary.csv"),
  row.names = FALSE
)

cat("Saved ground_truth_summary.{rds,csv}\n")

# -------------------------
# Quick sanity checks
# -------------------------

cat("\nSanity checks:\n")
cat("  Unique LR pairs:", nrow(GT_summary), "\n")
cat("  Max frequency:   ", max(GT_summary$freq), "\n")
cat("  Samples:         ", paste(names(lig_params_list), collapse = ", "), "\n")

cat("\nDone.\n")


#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
})

# ============================================================
# Config
# ============================================================

GT_PATH <- "simulations/ground_truth/ground_truth_summary.rds"
RUN_BASE_DIRS <- c(
  "simulations/tensor_runs",
  "simulations"
)

TOP_K <- 50
OUT_DIR <- "simulations/ground_truth"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load ground truth
# ============================================================

GT <- readRDS(GT_PATH) |>
  mutate(gt_key = paste(target, regulator, sep = "|"))

cat("Loaded ground truth with", nrow(GT), "unique LR pairs\n")

# ============================================================
# Helpers
# ============================================================

# ---- Extract genes from node name ----
# Format: 5_2_2_ITGA2_ITGB1  -> c("ITGA2","ITGB1")
extract_genes_from_node <- function(x) {
  parts <- strsplit(x, "_")[[1]]
  if (length(parts) < 4) return(character(0))
  gene_part <- paste(parts[4:length(parts)], collapse = "_")
  unlist(strsplit(gene_part, "_"))
}

# ---- Parse top_pos_interactions.csv ----
parse_top_pos <- function(csv_path) {
  df <- read.csv(csv_path, header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("source_node", "target_node", "score")
  
  parsed <- lapply(seq_len(nrow(df)), function(i) {
    ligs <- extract_genes_from_node(df$source_node[i])
    recs <- extract_genes_from_node(df$target_node[i])
    
    if (length(ligs) == 0 || length(recs) == 0) return(NULL)
    
    expand.grid(
      ligand = ligs,
      receptor = recs,
      score = df$score[i],
      rank = i,
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, parsed)
}

# ---- Evaluate top-K vs ground truth ----
evaluate_topk <- function(csv_path, GT, K = 50) {
  
  parsed <- parse_top_pos(csv_path)
  if (is.null(parsed) || nrow(parsed) == 0) return(NULL)
  
  parsed <- parsed |>
    mutate(key = paste(ligand, receptor, sep = "|"))
  
  topK <- parsed |>
    arrange(rank) |>
    filter(rank <= K)
  
  eval <- topK |>
    left_join(GT, by = c("key" = "gt_key"))
  
  data.frame(
    GT_hits_top50 = sum(!is.na(eval$freq)),
    GT_fraction_top50 = mean(!is.na(eval$freq)),
    Mean_GT_rank = mean(eval$rank[!is.na(eval$freq)], na.rm = TRUE),
    Mean_tensor_score_GT = mean(eval$score[!is.na(eval$freq)], na.rm = TRUE),
    Effect_weighted_recall_top50 =
      sum(eval$sum_effect[!is.na(eval$freq)], na.rm = TRUE) /
      sum(GT$sum_effect)
  )
}

# ============================================================
# Main loop
# ============================================================

results <- list()

for (base in RUN_BASE_DIRS) {
  if (!dir.exists(base)) next
  
  run_dirs <- list.dirs(base, recursive = FALSE)
  
  for (d in run_dirs) {
    csv <- file.path(d, "top_pos_interactions.csv")
    if (!file.exists(csv)) next
    
    cat("Evaluating:", d, "\n")
    
    stats <- evaluate_topk(csv, GT, TOP_K)
    if (is.null(stats)) next
    
    stats$run <- basename(d)
    stats$base <- base
    
    results[[length(results) + 1]] <- stats
  }
}

results_df <- bind_rows(results)

# ============================================================
# Save outputs
# ============================================================

write.csv(
  results_df,
  file.path(OUT_DIR, "top50_tensor_vs_ground_truth.csv"),
  row.names = FALSE
)

saveRDS(
  results_df,
  file.path(OUT_DIR, "top50_tensor_vs_ground_truth.rds")
)

cat("\nSaved results to:\n")
cat(" -", file.path(OUT_DIR, "top50_tensor_vs_ground_truth.csv"), "\n")
cat(" -", file.path(OUT_DIR, "top50_tensor_vs_ground_truth.rds"), "\n")
cat("\nDone.\n")
