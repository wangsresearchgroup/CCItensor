library(Seurat)
library(dplyr)
library(Azimuth)
library(CellChat)


# Input merged unfiltered scRNA-seq dataset and assign identities
merged<- readRDS("merged_seurat_object.rds")

n <- length(merged)
for(i in 1:n) {
  
  pdata <- merged[[i]]
  pdata[["percent.mt"]] <- PercentageFeatureSet(pdata, pattern = "^MT-")
  pdata <- subset(pdata, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  
  #Insert filtering here if not done already
  pdata <- RunAzimuth(pdata, reference = "lungref",k.weight = 40)
  p1 <- DimPlot(pdata, group.by = "predicted.ann_level_3", label = TRUE, label.size = 3) + NoLegend()
  #p2 <- DimPlot(pdata, group.by = "Method")
  p1
  
  pdata <- NormalizeData(pdata, normalization.method = "LogNormalize", scale.factor = 10000)
  pdata <- FindVariableFeatures(pdata, selection.method = "vst", nfeatures = 2000)
  pdata <- ScaleData(pdata)
  pdata <- RunPCA(pdata, features = VariableFeatures(object = pdata))
  pdata <- FindNeighbors(pdata, dims = 1:10)
  pdata <- FindClusters(pdata, resolution = 0.7)
  pdata <- RunUMAP(pdata, dims = 1:10)
  
  new.cluster.ids <- pdata$predicted.ann_level_3

  pdata$orig.ident <- pdata@active.ident

  p1 <- DimPlot(pdata, group.by = "predicted.ann_level_3", label = TRUE, label.size = 3) + NoLegend()
  plot1 <- DimPlot(pdata, reduction = "umap", pt.size = 0.5) + NoLegend()
  
  CellChat <- function(x) {
    
    
    result <- vector(mode = "list", length = 3)
    
    fdata <- x
    
    # CellChat Preparation
    cellchat <- createCellChat(object = fdata, group.by = "predicted.ann_level_3", assay = "RNA")
    
    
    CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
    showDatabaseCategory(CellChatDB)
    
    # Subset the expression data of signaling genes for saving computation cost
    cellchat@DB <- CellChatDB
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    #cellchat <- projectData(cellchat, PPI.human)
    
    # Compute the communication probability and infer cellular communication network
    ptm = Sys.time()
    cellchat <- computeCommunProb(cellchat, type = "triMean")
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    
    # Extract inferred cellular communication network as a data frame
    df.net <- subsetCommunication(cellchat)
    
    # Infer the cell-cell communication at signaling pathway level
    cellchat <- computeCommunProbPathway(cellchat)
    
    # Calculate the aggregated cell-cell communication network
    cellchat <- aggregateNet(cellchat)
    
    # Circle plot showing number of interactions or total interaction strength (weights) between any two cell groups
    ptm = Sys.time()
    groupSize <- as.numeric(table(cellchat@idents))
    par(mfrow = c(1, 2), xpd = TRUE)
    plot4 <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Number of interactions")
    plot5 <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Interaction weights/strength")
    mat <- cellchat@net$weight
    par(mfrow = c(3, 4), xpd = TRUE)
    result[[1]] <- df.net
    result[[2]] <- mat
    result[[3]] <- cellchat
    
    return(result)
  }
  
  res <- CellChat(pdata)
  merged[[i]] <- res
  print(i)
}

names(merged)<- identities

saveRDS(merged,file="your_output_cellchat_analysis")