ensure_packages <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message("Installing missing package: ", p)
      install.packages(
        p,
        repos = "https://cloud.r-project.org",
        lib = .libPaths()[1]
      )
    }
    library(p, character.only = TRUE)
  }
}

ensure_packages(c(
  "CellChat",
  "Seurat",
  "presto",
  "anndata",
  "SeuratObject",
  "tidyverse",
  "Matrix",
  "scMultiSim",
  "dplyr"
))

## ===============================================================
##  Load packages (assumed installed in rstudio-sc)
## ===============================================================
suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
  library(presto)
  library(anndata)
  library(SeuratObject)
  library(tidyverse)
  library(Matrix)
  library(scMultiSim)
  library(dplyr)
})

## ===============================================================
##  Load CellChat DB
## ===============================================================
data(CellChatDB.human)
db <- CellChatDB.human
rm(CellChatDB.human)

## ===============================================================
##  Utility: simulation parameter rules
## ===============================================================
get_sim_params <- function(sample_id) {
  if (sample_id <= 4) {
    list(numlink = 10, tree = Phyla3())
  } else if (sample_id <= 7) {
    list(numlink = 10, tree = Phyla5())
  } else {
    list(numlink = 20, tree = Phyla3())
  }
}

## ===============================================================
##  Generate ligand–receptor ground truth
## ===============================================================
generate_lig_params <- function(db, numlink, seed) {
  set.seed(seed)
  id <- sample(seq_len(nrow(db$interaction)), numlink)
  
  lig <- db$interaction$ligand.symbol[id]
  rec <- db$interaction$receptor.symbol[id]
  
  tar <- reg <- NULL
  for (i in seq_len(numlink)) {
    l1 <- unlist(strsplit(lig[i], ","))
    r1 <- unlist(strsplit(rec[i], ","))
    tar <- c(tar, rep(l1, each = length(r1)))
    reg <- c(reg, rep(r1, times = length(l1)))
  }
  
  eff <- runif(length(tar), min = 1, max = 10)
  
  data.frame(
    target = tar,
    regulator = reg,
    effect = eff
  )
}

## ===============================================================
##  Four-band wavelet (UNCHANGED)
## ===============================================================
four_band <- function(X) {
  n <- nrow(X); m <- ncol(X)
  h0 <- c(-0.067371764,0.094195111,0.40580489,0.567371764,
          0.567371764,0.40580489,0.094195111,-0.067371764)
  h1 <- c(-0.094195111,0.067371764,0.567371764,0.40580489,
          -0.40580489,-0.567371764,-0.067371764,0.094195111)
  h2 <- c(-0.094195111,-0.067371764,0.567371764,-0.40580489,
          -0.40580489,0.567371764,-0.067371764,-0.094195111)
  h3 <- c(-0.067371764,-0.094195111,0.40580489,-0.567371764,
          0.567371764,-0.40580489,0.094195111,0.067371764)
  
  nn <- 4 * ceiling(n / 4)
  vv1 <- ww1 <- ww2 <- ww3 <- matrix(0, nn, m)
  
  for (k in seq_len(m)) {
    ss <- X[, k]
    if (n %% 4) ss <- c(ss, rep(0, 4 - n %% 4))
    ss <- c(ss, ss[1:4])
    
    v1 <- w1 <- w2 <- w3 <- numeric(nn / 4)
    j <- 1
    for (i in seq_len(nn / 4)) {
      v1[i] <- sum(h0 * ss[j:(j+7)])
      w1[i] <- sum(h1 * ss[j:(j+7)])
      w2[i] <- sum(h2 * ss[j:(j+7)])
      w3[i] <- sum(h3 * ss[j:(j+7)])
      j <- j + 4
    }
    
    v1 <- c(v1[nn/4], v1)
    w1 <- c(w1[nn/4], w1)
    w2 <- c(w2[nn/4], w2)
    w3 <- c(w3[nn/4], w3)
    
    for (i in seq_len(nn / 4)) {
      idx <- (4*i-3):(4*i)
      vv1[idx,k] <- h0[5:8]*v1[i] + h0[1:4]*v1[i+1]
      ww1[idx,k] <- h1[5:8]*w1[i] + h1[1:4]*w1[i+1]
      ww2[idx,k] <- h2[5:8]*w2[i] + h2[1:4]*w2[i+1]
      ww3[idx,k] <- h3[5:8]*w3[i] + h3[1:4]*w3[i+1]
    }
  }
  
  list(
    low   = vv1[1:n,],
    high1 = ww1[1:n,],
    high2 = ww2[1:n,],
    high3 = ww3[1:n,]
  )
}

## ===============================================================
##  CellChat pipeline
## ===============================================================
cellchat_pipeline <- function(seurat_object, WT_matrix) {
  WT_matrix <- as.matrix(WT_matrix)
  rownames(WT_matrix) <- rownames(seurat_object)
  colnames(WT_matrix) <- colnames(seurat_object)
  
  input <- SetAssayData(seurat_object, "RNA", "data", WT_matrix)
  
  cellchat <- createCellChat(input, group.by = "cell.type", assay = "RNA")
  cellchat@DB <- db
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- aggregateNet(cellchat)
  
  list(
    df = subsetCommunication(cellchat),
    mat = cellchat@net$weight,
    cellchat = cellchat
  )
}

## ===============================================================
##  Run ALL simulations
## ===============================================================
BASE_OUT <- "simulations"
dir.create(BASE_OUT, showWarnings = FALSE)

for (s in 1:10) {
  
  cat("\n===== Running sample", s, "=====\n")
  params <- get_sim_params(s)
  
  outdir <- file.path(BASE_OUT, sprintf("sample_%02d", s))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  lig_params <- generate_lig_params(db, params$numlink, seed = 100 + s)
  write.csv(lig_params, file.path(outdir, "lig_params.csv"), row.names = FALSE)
  
  options <- list(
    rand.seed = s,
    GRN = NA,
    num.genes = 100,
    num.cells = 200,
    num.cifs = 50,
    cif.sigma = 0.5,
    tree = params$tree,
    diff.cif.fraction = 0.8,
    do.velocity = TRUE,
    cci = list(
      params = lig_params,
      max.neighbors = 4,
      cell.type.interaction = "random",
      step.size = 0.5
    )
  )
  
  sim <- sim_true_counts(options)
  
  seu <- CreateSeuratObject(sim$counts, meta.data = sim$cell_meta)
  seu <- NormalizeData(seu)
  
  data_mat <- as.matrix(GetAssayData(seu, layer = "data"))
  
  wavelet <- four_band(data_mat)
  wavelet <- lapply(wavelet, function(x) { x[x < 0] <- 0; x })
  
  results <- list(
    baseline = cellchat_pipeline(seu, data_mat),
    low      = cellchat_pipeline(seu, wavelet$low),
    high1    = cellchat_pipeline(seu, wavelet$high1),
    high2    = cellchat_pipeline(seu, wavelet$high2),
    high3    = cellchat_pipeline(seu, wavelet$high3)
  )
  
  for (nm in names(results)) {
    saveRDS(results[[nm]], file.path(outdir, paste0("cellchat_", nm, ".rds")))
  }
}
