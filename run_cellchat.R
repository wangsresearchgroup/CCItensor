args <- commandArgs(trailingOnly = TRUE)

p1_file  <- args[1]  # input .rds
out_file <- args[2]  # output .rds

# ---------- Load inputs ----------
p1 <- readRDS(p1_file)

p1_clean <- GetAssayData(
  p1,
  assay = "RNA",
  layer = "data"
)

# Apply wavelet and re-normalize
working_matrix <- as.matrix(p1_clean)

four_band <- function(X){
  n=dim(X)[1]
  m=dim(X)[2]
  h0<-(c(-0.067371764,0.094195111,0.40580489,0.567371764,0.567371764,0.40580489,0.094195111,-0.067371764))
  h1<-(c(-0.094195111,0.067371764, 0.567371764 ,0.40580489,-0.40580489,-0.567371764,-0.067371764,0.094195111))
  h2<-(c(-0.094195111,-0.067371764,0.567371764,-0.40580489,-0.40580489,0.567371764,-0.067371764,-0.094195111))
  h3<-(c(-0.067371764,-0.094195111,0.40580489,-0.567371764,0.567371764,-0.40580489,0.094195111,0.067371764))
  
  nn=4*ceiling(n/4)
  
  vv1=matrix(0,nn,m)
  ww1=matrix(0,nn,m)
  ww2=matrix(0,nn,m)
  ww3=matrix(0,nn,m)
  for (k in 1:m){
    ss=X[,k]
    if(n%%4){ss=c(ss,rep(0,4-n%%4))}
    ss=c(ss,ss[1:4])
    v1=rep(0,nn/4)
    w1=rep(0,nn/4)
    w2=rep(0,nn/4)
    w3=rep(0,nn/4)
    j=1
    for (i in 1:(nn/4)){
      v1[i]=sum(h0*ss[j:(j+7)])
      w1[i]=sum(h1*ss[j:(j+7)])
      w2[i]=sum(h2*ss[j:(j+7)])
      w3[i]=sum(h3*ss[j:(j+7)])
      j=j+4
    }
    v1=c(v1[nn/4],v1)
    w1=c(w1[nn/4],w1)
    w2=c(w2[nn/4],w2)
    w3=c(w3[nn/4],w3)
    for (i in 1:(nn/4)){
      vv1[4*i-3,k]=h0[5]*v1[i]+h0[1]*v1[i+1]
      vv1[4*i-2,k]=h0[6]*v1[i]+h0[2]*v1[i+1]
      vv1[4*i-1,k]=h0[7]*v1[i]+h0[3]*v1[i+1]
      vv1[4*i,k]=h0[8]*v1[i]+h0[4]*v1[i+1]
      ww1[4*i-3,k]=h1[5]*w1[i]+h1[1]*w1[i+1]
      ww1[4*i-2,k]=h1[6]*w1[i]+h1[2]*w1[i+1]
      ww1[4*i-1,k]=h1[7]*w1[i]+h1[3]*w1[i+1]
      ww1[4*i,k]=h1[8]*w1[i]+h1[4]*w1[i+1]
      ww2[4*i-3,k]=h2[5]*w2[i]+h2[1]*w2[i+1]
      ww2[4*i-2,k]=h2[6]*w2[i]+h2[2]*w2[i+1]
      ww2[4*i-1,k]=h2[7]*w2[i]+h2[3]*w2[i+1]
      ww2[4*i,k]=h2[8]*w2[i]+h2[4]*w2[i+1]
      ww3[4*i-3,k]=h3[5]*w3[i]+h3[1]*w3[i+1]
      ww3[4*i-2,k]=h3[6]*w3[i]+h3[2]*w3[i+1]
      ww3[4*i-1,k]=h3[7]*w3[i]+h3[3]*w3[i+1]
      ww3[4*i,k]=h3[8]*w3[i]+h3[4]*w3[i+1]
    }
  }
  vv1=vv1[1:n,]
  ww1=ww1[1:n,]
  ww2=ww2[1:n,]
  ww3=ww3[1:n,]
  
  list(low = vv1, high1 = ww1, high2 = ww2, high3 = ww3)
}

wavelet_result <- four_band(working_matrix)

wavelet_rescaled <- lapply(
  wavelet_result,
  function(mat) {
    mat <- as.matrix(mat)
    mat[mat < 0] <- 0
    mat
  }
)

cat("Wavelet finished\n")

# ---------- CellChat pipeline ----------
cellchat_pipeline <- function(seurat_object, WT_matrix) {
  WT_matrix <- as.matrix(WT_matrix)
  rownames(WT_matrix) <- rownames(seurat_object)
  colnames(WT_matrix) <- colnames(seurat_object)
  
  input <- SetAssayData(
    object = seurat_object,
    assay = "RNA",
    layer = "data",
    new.data = WT_matrix
  )
  
  cellchat <- createCellChat(
    object = input,
    group.by = "cell_type_predicted",
    assay = "RNA"
  )
  cellchat@DB <- CellChatDB.human
  
  cellchat <- subsetData(cellchat)
  
  if (nrow(cellchat@data.signaling) == 0) {
    stop("No signaling genes found.")
  }
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat@idents <- droplevels(cellchat@idents)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  
  df.net <- subsetCommunication(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  list(
    df = df.net,
    mat = cellchat@net$weight,
    obj = cellchat
  )
}
cellchat_pipeline<- function(seurat_object, WT_matrix) {
  # ---- Ensure matrix has correct dimnames ----
  WT_matrix <- as.matrix(WT_matrix)
  rownames(WT_matrix) <- rownames(seurat_object)
  colnames(WT_matrix) <- colnames(seurat_object)
  
  # ---- Meta must have rownames = cell barcodes ----
  meta <- seurat_object@meta.data
  meta <- meta[colnames(WT_matrix), , drop = FALSE]
  
  if (!"cell_type_predicted" %in% colnames(meta)) {
    stop("meta.data is missing column: cell_type_predicted")
  }
  
  # ---- Create CellChat from matrix+meta (NO Seurat GetAssayData(slot=) call) ----
  cellchat <- CellChat::createCellChat(
    object = WT_matrix,
    meta   = meta,
    group.by = "cell_type_predicted"
  )
  
  cellchat@DB <- CellChat::CellChatDB.human
  
  # critical: creates data.signaling
  cellchat <- CellChat::subsetData(cellchat)
  
  if (nrow(cellchat@data.signaling) == 0) stop("No signaling genes found.")
  
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat@idents <- droplevels(cellchat@idents)
  
  cellchat <- CellChat::computeCommunProb(cellchat, type = "triMean")
  df.net   <- CellChat::subsetCommunication(cellchat)
  
  cellchat <- CellChat::aggregateNet(cellchat)
  
  list(
    df  = df.net,
    mat = cellchat@net$weight,
    obj = cellchat
  )
}

# ---------- Run ----------
cat("Running untransformed\n")
p1_untransformed <- cellchat_pipeline(p1, working_matrix)

cat("Running low\n")
p1_approx <- cellchat_pipeline(p1, wavelet_rescaled$low)

cat("Running high1\n")
p1_h1 <- cellchat_pipeline(p1, wavelet_rescaled$high1)

cat("Running high2\n")
p1_h2 <- cellchat_pipeline(p1, wavelet_rescaled$high2)

cat("Running high3\n")
p1_h3 <- cellchat_pipeline(p1, wavelet_rescaled$high3)

# ---------- Save ----------
results <- list(
  approx = p1_approx,
  untransformed = p1_untransformed,
  h1 = p1_h1,
  h2 = p1_h2,
  h3 = p1_h3
)

saveRDS(results, out_file)

cat("CellChat pipeline finished successfully\n")
