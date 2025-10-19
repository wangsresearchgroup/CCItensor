library(reshape2)
library(ggplot2)
library(igraph)
# load metadata 
meta <- read.csv("CCI metadata - Sheet1.csv")

# metadata filtering
classification = ifelse(grepl("AS", meta$Sample), "squamous", "adenocarcinoma")
meta[,2] <- classification
colnames(meta)[2] <- "classification"
adeno_indices <- meta$classification == "adenocarcinoma"
squamous_indices <- meta$classification == "squamous"
adeno_meta <- meta[adeno_indices,]
names(adeno_meta)[1] <- "Sample"
squamous_meta <- meta[squamous_indices,]
names(squamous_meta)[1] <- "Sample"

merged_to_table_list <- function(merged_list) {
  # extract CCI tables from merged object in 
  l <- length(merged_list)
  for(i in 1:l) {
    merged_list[[i]] <- merged_list[[i]][[1]]
  }
  return(merged_list)
}

# total_table_list is a table of combined ligand-receptor tables with columns 'source', 'target', 'ligand', receptor', 'prob', 'pval', 'interaction_name', interaction_name_2', pathway_name', 'annotation'
total_table_list <- readRDS("merged_CCI_tables")


add_to_dict <- function(row, dict) {
  # function to check if an interaction described by the row(combination of source, target, ligand and receptor) matches an interaction in dict. If there is no match, the interaction is added to dict. 
  key <- paste(row[1:4], collapse = "_")
  if(key %in% names(dict)) {
    dict[[key]] <- dict[[key]] + 1
  }
  else {
    dict[[key]] <- 1
  }
  return(dict)
}
find_frequency <- function(table_list) {
  # Uses add_to_dict function to determine frequency of each interaction in the total table list
  frequency_dict <- list()
  for(table in table_list) {
    apply(table, MARGIN = 1, FUN = function(row) {
      frequency_dict <<- add_to_dict(row, frequency_dict)
    }
    )
  }
  return(frequency_dict)
}

freq <- find_frequency(total_table_list)

filter_tables <- function (table_list, frequency_dict, threshold) {
  # removes interactions that appear less than threshold frequency in entire sample
  filtered_frequency <- frequency_dict[sapply(frequency_dict, function(x) {x >= threshold})]
  
  filtered_table_list <- list()
  for (i in seq_along(table_list)) {
    table <- table_list[[i]]
    table <- table[apply(table, MARGIN = 1, function(row) {
      paste(row[1:4], collapse = "_") %in% names(filtered_frequency)
    }), ]
    filtered_table_list[[i]] <- table
  }
  names(filtered_table_list) <- names(table_list)
  return(filtered_table_list)
}

filtered_table_list <- filter_tables(total_table_list, find_frequency(total_table_list), 45)
filtered_adeno_list <- filtered_table_list[adeno_indices]
filtered_squamous_list <- filtered_table_list[squamous_indices]


five_tensor <- function(input_list) {
  # transforms CCI table list into a five-way tensor with source, target, ligand, receptor and sample axes
  combined_table <- as.data.frame(input_list[[1]])
  for(table in input_list) {
    combined_table <- rbind(combined_table, table)
  }
  sources <- unique(combined_table$source)
  targets <- unique(combined_table$target)
  ligand <- unique(combined_table$ligand)
  receptor <- unique(combined_table$receptor)
  
  tensor <- array(0, dim = c(length(sources), length(targets), length(ligand), length(receptor), length(input_list)),
                  dimnames = list(sources, targets, ligand, receptor, names(input_list)))
  
  for(table_idx in seq_along(input_list)) {
    table <- input_list[[table_idx]]
    samp <- names(input_list)[table_idx]
    for (i in seq_len(nrow(table))) {
      src <- as.character(table$source[i])
      tgt <- as.character(table$target[i])
      lig <- as.character(table$ligand[i])
      rec <- as.character(table$receptor[i])
      
      prob <- table$prob[i]
      
      # Fill the tensor at the appropriate index
      tensor[src, tgt, lig, rec, samp] <- prob
    }
    
    
  }
  return(list(tensor, sources, targets, ligand, receptor, sample))
}
tensor_res_A <- five_tensor(filtered_adeno_list)
tensor_res_S <- five_tensor(filtered_squamous_list)

library(rTensor)
# tensor decomposition
tensor_5_A <- as.tensor(tensor_res_A[[1]])
tensor_5_S <- as.tensor(tensor_res_S[[1]])

# CP tensor decomposition for LUAD and LUSC tensors
ntf_res_A <- cp(tensor_5_A, 80, 50)
ntf_res_S <- cp(tensor_5_S, 80, 50)

saveRDS(ntf_res_A, "t_decomp_80,50_a")
saveRDS(ntf_res_S, "t_decomp_80,50_s")

ntf_res_A  <- readRDS("t_decomp_80,50_a")
ntf_res_S <- readRDS("t_decomp_80,50_s")

# function form 
tensor_decomp <- function(filtered_list, num_components, num_iter) {
  tensor_res <- five_tensor(filtered_list)
  tensor_5 <- as.tensor(tensor_res[[1]])
  decomp_result <- cp(tensor_5, num_components, num_iter)
  return(decomp_result, tensor_res)
}

# Compute the outer product
compute_outer <- function(vec1, vec2) {
  # computes outer product of vectors from factor matrices
  outer_product <- array(
    outer(vec1, vec2),
    dim = c(length(vec1), length(vec2)),
    dimnames = list(names(vec1), names(vec2))
  )
  return(outer_product)
}

construct_two_array <- function(mat_1, mat_2) {
  # create list of factor matrices corresponding to the otuer product of mat1 and mat2 factors
  # Initialize an empty list with the required length
  array_list <- vector("list", ncol(mat_1))
  
  # Loop through columns
  for (i in seq_len(ncol(mat_1))) {
    # Compute the outer product for each set of vectors
    array_list[[i]] <- compute_outer(mat_1[,i], mat_2[,i])
  }
  
  return(array_list)
}

flatten_array <- function(array) {
  # flattens array into a vector 
  
  # Get all combinations of dimension names
  index_combinations <- expand.grid(dimnames(array))
  
  # Create names by concatenating dimension names with underscores
  names_flat <- apply(index_combinations, 1, paste, collapse = "_")
  
  # Flatten the array into a vector
  flat_array <- as.vector(array)
  
  # Assign names to the flattened vector
  names(flat_array) <- names_flat
  
  return(flat_array)
}



create_flat_list <- function(array_list) {
  # use flatten_array function to turn list of arrays into flattened vectors
  flat_list <- vector("list", length(array_list))
  for(i in seq_len(length(array_list))) {
    flat_list[[i]] <- flatten_array(array_list[[i]])
  }
  return(flat_list)
}


flat_list_to_flat_table <- function(flat_list) {
  # combine vectors in flat list into a table
  table <- do.call(cbind, flat_list)
  return(table)
}

compute_inner <- function(s_l_table, t_r_table) {
  # compute inner product  matrix given sample matrix and combined flat matrix
  inner_matrix <- matrix(0, nrow(s_l_table),nrow(t_r_table))
  for(i in seq_len(nrow(s_l_table))) {
    for(j in seq_len(nrow(t_r_table))) {
      inner_matrix[i, j] <- sum(s_l_table[i, ] * t_r_table[j, ])
      
    }
  }
  rownames(inner_matrix) <- rownames(s_l_table)
  colnames(inner_matrix) <- rownames(t_r_table)
  inner_matrix <- t(inner_matrix)
  return(inner_matrix)
}


four_band <- function(x){
  # computes four-band DWT on matrix x
  l <- length(x)
  fourbandlow<- vector("list",l)
  fourbandhigh1<- vector("list",l)
  fourbandhigh2<- vector("list",l)
  fourbandhigh3<- vector("list",l)
  for(p in 1:l){
    X <- x[[p]]
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
    
    fourbandlow[[p]] <- vv1
    fourbandhigh1[[p]]<-ww1
    fourbandhigh2[[p]]<-ww2
    fourbandhigh3[[p]]<-ww3
  }
  result= list(fourbandlow,fourbandhigh1,fourbandhigh2,fourbandhigh3)
  return (result)
}
array_list_to_flat_table <- function(array_list) {
  # converts list of arrays to one flat table where each array is flattened to a column
  return(flat_list_to_flat_table(create_flat_list(array_list)))
}
inner_pipeline <- function(ntf_res, tensor_res, sample_names) {
  # generates interaction matrix by taking outer product of source + ligand factor matrices and target + receptor matrices,
  # and taking the inner product of each source_ligand + target_receptor pair 
  source_matrix <- ntf_res$U[[1]] 
  rownames(source_matrix) <- tensor_res[[2]]
  target_matrix <- ntf_res$U[[2]]
  rownames(target_matrix) <- tensor_res[[3]]
  ligand_matrix <- ntf_res$U[[3]]
  rownames(ligand_matrix) <- tensor_res[[4]]
  receptor_matrix <- ntf_res$U[[4]]  
  rownames(receptor_matrix) <- tensor_res[[5]]
  sample_matrix <- ntf_res$U[[5]]
  rownames(sample_matrix) <- sample_names
  
  s_l_array_list <- construct_two_array(source_matrix, ligand_matrix)
  t_r_array_list <- construct_two_array(target_matrix, receptor_matrix)

  s_l_four_band_result <- four_band(s_l_array_list)
  s_l_array_list_approx <- s_l_four_band_result[[1]]
  s_l_array_list_h1 <- s_l_four_band_result[[2]]
  s_l_array_list_h2 <- s_l_four_band_result[[3]]
  s_l_array_list_h3 <- s_l_four_band_result[[4]]
  
  t_r_four_band_result <- four_band(t_r_array_list)
  t_r_array_list_approx <- t_r_four_band_result[[1]]
  t_r_array_list_h1 <- t_r_four_band_result[[2]]
  t_r_array_list_h2 <- t_r_four_band_result[[3]]
  t_r_array_list_h3 <- t_r_four_band_result[[4]]
  
  s_l_table <- array_list_to_flat_table(s_l_array_list)
  t_r_table <- array_list_to_flat_table(t_r_array_list)
  s_l_table_approx <- array_list_to_flat_table(s_l_array_list_approx)
  s_l_table_h1 <- array_list_to_flat_table(s_l_array_list_h1)
  s_l_table_h2 <- array_list_to_flat_table(s_l_array_list_h2)
  s_l_table_h3 <- array_list_to_flat_table(s_l_array_list_h3)
  t_r_table_approx <- array_list_to_flat_table(t_r_array_list_approx)
  t_r_table_h1 <- array_list_to_flat_table(t_r_array_list_h1)
  t_r_table_h2 <- array_list_to_flat_table(t_r_array_list_h2)
  t_r_table_h3 <- array_list_to_flat_table(t_r_array_list_h3)
  
  inner_matrix <- compute_inner(s_l_table, t_r_table)
  print(inner_matrix)
  inner_matrix_approx <- compute_inner(s_l_table_approx, t_r_table_approx)
  inner_matrix_h1 <- compute_inner(s_l_table_h1, t_r_table_h1)
  inner_matrix_h2 <- compute_inner(s_l_table_h2, t_r_table_h2)
  inner_matrix_h3 <- compute_inner(s_l_table_h3, t_r_table_h3)
  
  rownames(inner_matrix_approx) <- rownames(inner_matrix)
  colnames(inner_matrix_approx) <- colnames(inner_matrix)
  rownames(inner_matrix_h1) <- rownames(inner_matrix)
  colnames(inner_matrix_h1) <- colnames(inner_matrix)
  rownames(inner_matrix_h2) <- rownames(inner_matrix)
  colnames(inner_matrix_h2) <- colnames(inner_matrix)
  rownames(inner_matrix_h3) <- rownames(inner_matrix)
  colnames(inner_matrix_h3) <- colnames(inner_matrix)
  return(list(
    inner = inner_matrix,
    approx_inner = inner_matrix_approx,
    h1_inner = inner_matrix_h1,
    h2_inner = inner_matrix_h2,
    h3_inner = inner_matrix_h3
  ))
}
 
# compute inner product interaction matrices
inner_matrix_A <- inner_pipeline(ntf_res_A, tensor_res_A, adeno_meta$Sample)[[1]]
inner_matrix_A_approx <- inner_pipeline(ntf_res_A, tensor_res_A, adeno_meta$Sample)[[2]]
inner_matrix_A_h1 <- inner_pipeline(ntf_res_A, tensor_res_A, adeno_meta$Sample)[[3]]
inner_matrix_A_h2<- inner_pipeline(ntf_res_A, tensor_res_A, adeno_meta$Sample)[[4]]
inner_matrix_A_h3 <- inner_pipeline(ntf_res_A, tensor_res_A, adeno_meta$Sample)[[5]]

inner_matrix_S <- inner_pipeline(ntf_res_S, tensor_res_S, squamous_meta$Sample)[[1]]
inner_matrix_S_approx <- inner_pipeline(ntf_res_S, tensor_res_S, squamous_meta$Sample)[[2]]
inner_matrix_S_h1 <- inner_pipeline(ntf_res_S, tensor_res_S, squamous_meta$Sample)[[3]]
inner_matrix_S_h2 <- inner_pipeline(ntf_res_S, tensor_res_S, squamous_meta$Sample)[[4]]
inner_matrix_S_h3 <- inner_pipeline(ntf_res_S, tensor_res_S, squamous_meta$Sample)[[5]]

# Find percent zeros
calculate_zero_percentage <- function(matrix_input) {
  # Count the total number of elements in the matrix
  matrix_input <- as.matrix(matrix_input)
  total_elements <- length(matrix_input)
  
  # Count the number of zero elements in the matrix
  zero_elements <- sum(abs(matrix_input) < 10^(-10))
  
  # Calculate the percentage of zero elements
  zero_percentage <- (zero_elements / total_elements) * 100
  
  # Return the percentage of zero elements
  return(zero_percentage)
}
percent_zero_inner_A <- calculate_zero_percentage(inner_matrix_A)
percent_zero_inner_A_approx <- calculate_zero_percentage(inner_matrix_A_approx)
percent_zero_inner_A_h1 <- calculate_zero_percentage(inner_matrix_A_h1)
percent_zero_inner_A_h2 <- calculate_zero_percentage(inner_matrix_A_h2)
percent_zero_inner_A_h3 <- calculate_zero_percentage(inner_matrix_A_h3)



plot_inner_heatmap <- function(inner_matrix) {
  # plot heatmaps of inner product matrices 
  mat_melted <- melt(inner_matrix, varnames = c("Row", "Column"))
  
  plot <- ggplot(mat_melted, aes(x = Column, y = Row, fill = value)) +
    geom_tile() +

    scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, 
    ) + 
    theme_minimal() + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 5))
  
  return(plot)
}

add_names <- function(w_matrix, o_matrix) {
  rownames(w_matrix) <- rownames(o_matrix)
  colnames(w_matrix) <- colnames(o_matrix)
  return(w_matrix)
}

plotA <- plot_inner_heatmap(inner_matrix_A)
plotA_approx <- plot_inner_heatmap(add_names(inner_matrix_A_approx, inner_matrix_A))
plotA_h1 <- plot_inner_heatmap(add_names(inner_matrix_A_h1, inner_matrix_A))
plotA_h2 <- plot_inner_heatmap(add_names(inner_matrix_A_h2, inner_matrix_A))
plotA_h3 <- plot_inner_heatmap(add_names(inner_matrix_A_h3, inner_matrix_A))

plotS <- plot_inner_heatmap(inner_matrix_S)
plotS_approx <- plot_inner_heatmap(add_names(inner_matrix_S_approx, inner_matrix_S))
plotS_h1 <- plot_inner_heatmap(add_names(inner_matrix_S_h1, inner_matrix_S))
plotS_h2 <- plot_inner_heatmap(add_names(inner_matrix_S_h2, inner_matrix_S))
plotS_h3 <- plot_inner_heatmap(add_names(inner_matrix_S_h3, inner_matrix_S))

png("inner_heatmap_A.png",width=2000,height=1500,res=600)
print(plotA)
dev.off()

png("inner_heatmap_S.png",width=1000,height=1000,res=300)
print(plotS)
dev.off()

png("inner_heatmap_A_approx.png",width=2000,height=1500,res=600)

print(plotA_approx)
dev.off()

png("inner_heatmap_A_h1.png",width=2000,height=1500,res=600)
print(plotA_h1)
dev.off()

png("inner_heatmap_A_h2.png",width=2000,height=1500,res=600)
print(plotA_h2)
dev.off()

png("inner_heatmap_A_h3.png",width=2000,height=1500,res=600)
print(plotA_h3)
dev.off()

png("inner_heatmap_S_approx.png",width=1000,height=1000,res=300)
print(plotS_approx)
dev.off()

png("inner_heatmap_S_h1.png",width=1000,height=1000,res=300)
print(plotS_h1)
dev.off()

png("inner_heatmap_S_h2.png",width=1000,height=1000,res=300)
print(plotS_h2)
dev.off()

png("inner_heatmap_S_h3.png",width=1000,height=1000,res=300)
print(plotS_h3)
dev.off()

# graph generation
create_adjacency_matrix <- function(mat, threshold) {
  # creates adjacency matrix of top interactions above threshold percentge
  maximum <- max(as.vector(mat))
  abs_thresh <- maximum * threshold
  for(i in 1:nrow(mat)) {
    for(j in 1:ncol(mat)) {
      if(mat[i,j] >= abs_thresh ) {
        mat[i,j] <- 1
      } else {
        mat[i,j] <- 0
      }
    }
  }
  return(mat)
}
# define color palette for cell types in network graph
vertex_cell_types <- lapply(c(colnames(inner_matrix_A),rownames(inner_matrix_A)), 
                            function(x) gsub("_.*$", "", x))

vertex_l_r <- lapply(c(colnames(inner_matrix_A),rownames(inner_matrix_A)), 
                     function(x) gsub("$_", "", x))

#vertex_colors <- brewer.pal(length(unique(vertex_cell_types)), "Dark2")
vertex_colors <- c(
  "yellowgreen",  # Dark blue
  "#67153d",  # Deep 
  "darkgreen",  # Dark green
  "#FF4500",  # Dark orange red
  "#8B0000",  # Dark red
  "#4682B4",  # Steel blue
  "#9400D3",  # Purple
  "#6d4e3e",  # brown
  "darkorange",  # Tomato red
  "pink",  # Rich purple
  "black"   # Deep plum
)

names(vertex_colors) <- unique(vertex_cell_types)


create_graph <- function(inner_matrix, vertex_colors) {
  # create network graph using inner matrix and vertex colors 
  adj <- create_adjacency_matrix(inner_matrix, 0.1)
  u_adj <-create_adjacency_matrix(inner_matrix, 0)
  graph <- graph_from_biadjacency_matrix(adj, mode = "in", directed = TRUE, weighted = TRUE)
  u_graph <- graph_from_biadjacency_matrix(inner_matrix, mode = "in", directed = TRUE, weighted = TRUE)
  # Remove vertices with no edges
  graph.attributes <- vertex_attr_names(graph)
  vertex_cell_types <- lapply(V(graph)$name, 
                              function(x) gsub("_.*$", "", x))
  l_r_names <- lapply(V(graph)$name, 
                              function(x) gsub("^[^_]*_", "", x))
  # Apply vertex colors
  color_list <- sapply(vertex_cell_types, function(x) {
    for(i in 1:length(vertex_colors)) {
      index <- match(x, names(vertex_colors))
      return(vertex_colors[[index]])
    }
  })
  
  edge_list <- get.data.frame(graph, what = "edges")
  edge_weights <- sapply(1:nrow(edge_list), function(i) {
    weight <- inner_matrix[edge_list[i, 2], edge_list[i, 1]]
    return(weight) 
  })
  
  edge_list_u<- get.data.frame(u_graph, what = "edges")
  
  E(graph)$weight <- edge_weights
  V(graph)$size <- 2
  V(graph)$label <- l_r_names
  V(graph)$label.cex <- 1.2
  V(graph)$label.color <- color_list
  V(graph)$color <- "lightyellow"
  
  centrality <- eigen_centrality(u_graph)$vector
  V(graph)$size <- ifelse(centrality > 0, centrality * 40, 0)
  graph <- delete_vertices(graph, degree(graph) == 0)
  
  
 
  return(list(graph, edge_list,l_r_names))
  
}
graph_A <- create_graph(inner_matrix_A, vertex_colors)
graph_S <- create_graph(inner_matrix_S, vertex_colors)

graph_A_approx <- create_graph(inner_matrix_A_approx, vertex_colors)
graph_S_approx <- create_graph(inner_matrix_S_approx, vertex_colors)

graph_A_h1 <- create_graph(inner_matrix_A_h1, vertex_colors)
graph_S_h1 <- create_graph(inner_matrix_S_h1, vertex_colors)

graph_A_h2 <- create_graph(inner_matrix_A_h2, vertex_colors)
graph_S_h2 <- create_graph(inner_matrix_S_h2, vertex_colors)

graph_A_h3 <- create_graph(inner_matrix_A_h3, vertex_colors)
graph_S_h3 <- create_graph(inner_matrix_S_h3, vertex_colors)

png("inner_graph_A.png",width=3100,height=3100,res=300)
plot(graph_A[[1]])
dev.off()

png("inner_graph_S.png",width=3000,height=3000,res=300)
plot(graph_S[[1]])
dev.off()

png("inner_graph_A_approx.png",width=3100,height=3100,res=300)
plot(graph_A_approx[[1]])
dev.off()

png("inner_graph_S_approx.png",width=3000,height=3000,res=300)
plot(graph_S_approx[[1]])
dev.off()

png("inner_graph_A_h1.png",width=3100,height=3100,res=300)
plot(graph_A_h1[[1]])
dev.off()

png("inner_graph_S_h1.png",width=3000,height=3000,res=300)
plot(graph_S_h1[[1]])
dev.off()

png("inner_graph_A_h2.png",width=3100,height=3100,res=300)
plot(graph_A_h2[[1]])
dev.off()

png("inner_graph_S_h2.png",width=3000,height=3000,res=300)
plot(graph_S_h2[[1]])
dev.off()

png("inner_graph_A_h3.png",width=3100,height=3100,res=300)
plot(graph_A_h3[[1]])
dev.off()

png("inner_graph_S_h3.png",width=3000,height=3000,res=300)
plot(graph_S_h3[[1]])
dev.off()

# display legend
legend("topright", 
       legend = names(vertex_colors), 
       col = vertex_colors, 
       pch = c(19), 
       title = "Legend")


png("legend_graph_network.png",width=2500,height=2500,res=300)
plot(graph_A[[1]])
legend("topright", 
       legend = names(vertex_colors), 
       col = vertex_colors, 
       pch = c(19), 
       title = "Legend")
dev.off()

# identify top interactions
generate_top_signatures <- function(inner_matrix) {
  # generates top interaction with highest values in interaction matrix 
  inner_matrix[is.na(inner_matrix)] <- 0
  long_inner <- data.frame(
    ligand = rep(colnames(inner_matrix), each = nrow(inner_matrix)), 
    receptor = rep(rownames(inner_matrix), times = ncol(inner_matrix)), 
    interaction_value = as.vector(inner_matrix) 
  )
  
  ordered_inner <- long_inner[order(long_inner$interaction_value),] 
  threshold <- ceiling(nrow(ordered_inner)*0.0002)
  top_pos_interactions <- ordered_inner[(nrow(ordered_inner) - threshold + 1):nrow(ordered_inner),]
  top_pos_interactions <- top_pos_interactions[nrow(top_pos_interactions):1, , drop = FALSE]
  return(top_pos_interactions)
}

A_sigs <- generate_top_signatures(inner_matrix_A)
S_sigs <- generate_top_signatures(inner_matrix_S)
# Wavelet signatures generated from average of wavelet matrices

wav_A_average <- (inner_matrix_A_approx + inner_matrix_A_h1 + inner_matrix_A_h2 + inner_matrix_A_h3) / 4
rownames(wav_A_average) <- rownames(inner_matrix_A)
colnames(wav_A_average) <- colnames(inner_matrix_A)
wav_A_sigs <- generate_top_signatures(wav_A_average)

wav_S_average <- (inner_matrix_S_approx + inner_matrix_S_h1 + inner_matrix_S_h2 + inner_matrix_S_h3) / 4
rownames(wav_S_average) <- rownames(inner_matrix_S)
colnames(wav_S_average) <- colnames(inner_matrix_S)
wav_S_sigs <- generate_top_signatures(wav_S_average)

write.csv(A_sigs,"A_top_pos_interactions.csv", row.names = FALSE)
write.csv(wav_A_sigs,"A_top_wav_interactions.csv", row.names = FALSE)
write.csv(S_sigs,"S_top_pos_interactions.csv", row.names = FALSE)
write.csv(wav_S_sigs,"S_top_wav_interactions.csv", row.names = FALSE)

# umap based on tensor decomposoed sample factors 
UMAP <- function(factors, meta) {
  set.seed(100)
  umap_res<- umap(factors,n_neighbors=7)
  
  umap_data<- umap_res$layout
  
  # Perform k-means clustering on the UMAP embeddings
  set.seed(100)
  num_clusters <- 3# Adjust the number of clusters as needed
  umap_clusters<- kmeans(umap_data, centers = num_clusters)$cluster
  
  age <- meta$Age
  gender <- meta$Gender
  smoking <- meta$Smoking
  stage <- meta$Stages
  
  plot_umap <- function(umap_df, fac, color_pal) {
    ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = fac)) +
      geom_point(size = 3, alpha = 0.7) +  # Set point color to white
      theme_classic() +
      theme(legend.position = "right") + 
      scale_color_manual(values =  color_pal) +
      stat_ellipse(aes(group = umap_clusters), color = "black", linetype = "dashed") +  # First set of ellipses
      stat_ellipse(aes(color = fac))  # Second set of ellipses
    
  }
  
  umap_df <- data.frame(UMAP1 = umap_data[, 1], UMAP2 = umap_data[, 2], Cluster = as.factor(umap_clusters))
  
  plot_cluster <- plot_umap(umap_df, as.factor(umap_clusters), c("#F8766D","#00BA38", "#619CFF"))
  plot_gender <- plot_umap(umap_df, as.factor(gender), c("#F8766D", "#00BFC4"))
  plot_stage <- plot_umap(umap_df, as.factor(stage), c("green","purple"))
  
  
  return(list(plot_cluster,umap_clusters, plot_gender, plot_stage))
}
dimreduce <- function(tensor_decomp_res,meta) {
  factors_mode_5_u <- tensor_decomp_res$U[[5]]
  umap_mode5 <- UMAP(factors_mode_5_u, meta)
  
  return(umap_mode5)
}
library(ggplot2)
library(ggrepel)
library(umap)
umap_A <- dimreduce(ntf_res_A, meta[adeno_indices,])
umap_S <- dimreduce(ntf_res_S, meta[squamous_indices,])


png("umap_A.png",width=1500,height=1500,res=300)
umap_A[[1]]
dev.off()

png("umap_A_gender.png",width=1500,height=1500,res=300)
umap_A[[4]]
dev.off()

png("umap_A_stage.png",width=1500,height=1500,res=300)
umap_A[[5]]
dev.off()

png("umap_S_png",width=1500,height=1500,res=300)
umap_S[[1]]
dev.off()

# Fisher's exact test for umap clusters
adeno_meta <- meta[adeno_indices,]
cluster1 <- meta[adeno_indices,][umap_A[[2]] == 1,]
cluster2 <- meta[adeno_indices,][umap_A[[2]] == 2,]
cluster3 <- meta[adeno_indices,][umap_A[[2]] == 3,]

fisher_test <- function(label1, label2, category) {
  # generate fisher's exact test p_value for distribution of label1 and label2 amoong LUAD clusters
  total_l1 <- nrow(adeno_meta[adeno_meta[[category]] == label1,])
  total_l2 <- nrow(adeno_meta[adeno_meta[[category]] == label2,])
  
  cluster1_l1 <- nrow(cluster1[cluster1[[category]] == label1,])
  cluster1_l2 <- nrow(cluster1[cluster1[[category]] == label2,])
  cluster2_l1<- nrow(cluster2[cluster2[[category]] == label1,])
  cluster2_l2 <- nrow(cluster2[cluster2[[category]] == label2,])
  cluster3_l1 <- nrow(cluster3[cluster3[[category]] == label1,])
  cluster3_l2 <- nrow(cluster3[cluster3[[category]] == label2,])
  
  contigency_table <- rbind(c(cluster1_l1, cluster1_l2),
                            c(cluster2_l1, cluster2_l2),
                            c(cluster3_l1, cluster3_l2),
                            c(total_l1, total_l2)
  )
  rownames(contigency_table) <- c("1","2","3","total")
  colnames(contigency_table) <- c(label1, label2)
  fisher_input<- contigency_table[1:3,]
  fisher_result <- fisher.test(fisher_input)
  return(fisher_result)
}
fisher_result_stage <- fisher_test("I/II","III/IV","Stages")$p.value
fisher_result_gender <- fisher_test("M","F","Gender")$p.value

# subtype analysis 
a_meta <- meta[adeno_indices,]
a_meta$cluster <- umap_A[[2]]
a_c1_indices <- a_meta$cluster == 1
a_c2_indices <- a_meta$cluster == 2
a_c3_indices <- a_meta$cluster == 3
a_c1_names <- a_meta[a_c1_indices,]$Sample
a_c2_names <- a_meta[a_c2_indices,]$Sample
a_c3_names <- a_meta[a_c3_indices,]$Sample

s_meta <- meta[squamous_indices,]
s_meta$cluster <- umap_S[[2]]
s_c1_indices <- s_meta$cluster == 1
s_c2_indices <- s_meta$cluster == 2
s_c3_indices <- s_meta$cluster == 3
s_c1_names <- s_meta[s_c1_indices,]$Sample
s_c2_names <- s_meta[s_c2_indices,]$Sample
s_c3_names <- s_meta[s_c3_indices,]$Sample

filtered_list_a_c1 <- filtered_adeno_list[a_c1_indices]
filtered_list_a_c2 <- filtered_adeno_list[a_c2_indices]
filtered_list_a_c3 <- filtered_adeno_list[a_c3_indices]

filtered_list_s_c1 <- filtered_squamous_list[s_c1_indices]
filtered_list_s_c2 <- filtered_squamous_list[s_c2_indices]
filtered_list_s_c3 <- filtered_squamous_list[s_c3_indices]

tensor_res_A_c1 <- five_tensor(filtered_list_a_c1)
tensor_res_A_c2 <- five_tensor(filtered_list_a_c2)
tensor_res_A_c3 <- five_tensor(filtered_list_a_c3)

tensor_res_S_c1 <- five_tensor(filtered_list_s_c1)
tensor_res_S_c2 <- five_tensor(filtered_list_s_c2)
tensor_res_S_c3 <- five_tensor(filtered_list_s_c3)

# tensor decomposition
tensor_5_A_c1 <- as.tensor(tensor_res_A_c1[[1]])
tensor_5_A_c2 <- as.tensor(tensor_res_A_c2[[1]])
tensor_5_A_c3 <- as.tensor(tensor_res_A_c3[[1]])

tensor_5_S_c1 <- as.tensor(tensor_res_S_c1[[1]])
tensor_5_S_c3 <- as.tensor(tensor_res_S_c3[[1]])
tensor_5_S_c3 <- as.tensor(tensor_res_S_c3[[1]])



ntf_res_A_c1 <- tensor_decomp(filtered_list_a_c1, 80, 50)

ntf_res_A_c2 <- tensor_decomp(filtered_list_a_c2, 80, 50)

ntf_res_A_c3 <- tensor_decomp(filtered_list_a_c3, 30, 50)

ntf_res_S_c1 <- tensor_decomp(filtered_list_s_c1, 30, 50)
ntf_res_S_c2 <- tensor_decomp(filtered_list_s_c2, 80, 50)
ntf_res_S_c3 <- tensor_decomp(filtered_list_s_c3, 30, 50)

saveRDS(ntf_res_A_c1, "tensor_decomp_a_c1")
saveRDS(ntf_res_A_c2, "tensor_decomp_a_c2")
saveRDS(ntf_res_A_c3, "tensor_decomp_a_c3")
saveRDS(ntf_res_S_c1, "tensor_decomp_s_c1")
saveRDS(ntf_res_S_c2, "tensor_decomp_s_c2")
saveRDS(ntf_res_S_c3, "tensor_decomp_s_c3")

ntf_res_A_c1 <- readRDS("tensor_decomp_a_c1")
ntf_res_A_c2 <- readRDS("tensor_decomp_a_c2")
ntf_res_A_c3 <- readRDS("tensor_decomp_a_c3")

inner_results_a_c1 <- inner_pipeline(ntf_res_A_c1, tensor_res_A_c1, a_c1_names)
inner_results_a_c2 <- inner_pipeline(ntf_res_A_c2, tensor_res_A_c2, a_c2_names)
inner_results_a_c3 <- inner_pipeline(ntf_res_A_c3, tensor_res_A_c3, a_c3_names)
inner_results_s_c1 <- inner_pipeline(ntf_res_S_c1, tensor_res_S_c1, s_c1_names)
inner_results_s_c2 <- inner_pipeline(ntf_res_S_c2, tensor_res_S_c2, s_c2_names)
inner_results_s_c3 <- inner_pipeline(ntf_res_S_c3, tensor_res_S_c3, s_c3_names)
inner_A_c1 <- inner_results_a_c1[[1]]
inner_A_c2 <- inner_results_a_c2[[1]]
inner_A_c3 <- inner_results_a_c3[[1]]

inner_c1_S <- inner_results_s_c1[[1]]
inner_c2_S <- inner_results_s_c2[[1]]
inner_c3_S <- inner_results_s_c3[[1]]

inner_A_c1_approx <- inner_results_a_c1[[2]]
inner_A_c2_approx<- inner_results_a_c2[[2]]
inner_A_c3_approx <- inner_results_a_c3[[2]]

inner_A_c1_h1 <- inner_results_a_c1[[3]]
inner_A_c2_h1<- inner_results_a_c2[[3]]
inner_A_c3_h1 <- inner_results_a_c3[[3]]

inner_A_c1_h2 <- inner_results_a_c1[[4]]
inner_A_c2_h2<- inner_results_a_c2[[4]]
inner_A_c3_h2 <- inner_results_a_c3[[4]]

inner_A_c1_h3 <- inner_results_a_c1[[5]]
inner_A_c2_h3<- inner_results_a_c2[[5]]
inner_A_c3_h3 <- inner_results_a_c3[[5]]

custom_palette <- colorRampPalette(c("blue", "white", "red"))(20)

## subtype heatmaps

plotA_c1 <- plot_inner_heatmap(inner_A_c1)
plotA_c2 <- plot_inner_heatmap(inner_A_c2)
plotA_c3 <- plot_inner_heatmap(inner_A_c3)

plotS_c1 <- plot_inner_heatmap(inner_c1_S)
plotS_c2 <- plot_inner_heatmap(inner_c2_S)
plotS_c3 <- plot_inner_heatmap(inner_c3_S)


png("inner_heatmap_A_c1.png",width=1000,height=1000,res=300)
print(plotA_c1)
dev.off()

png("inner_heatmap_A_c2.png",width=1000,height=1000,res=300)
print(plotA_c2)
dev.off()

png("inner_heatmap_A_c3.png",width=1000,height=1000,res=300)
print(plotA_c3)
dev.off()

png("inner_heatmap_S_c1.png",width=1000,height=1000,res=300)
print(plotS_c1)
dev.off()

png("inner_heatmap_S_c2.png",width=1000,height=1000,res=300)
print(plotS_c2)
dev.off()

png("inner_heatmap_S_c3.png",width=1000,height=1000,res=300)
print(plotS_c3)
dev.off()


# Subtype network graphs
graph_A_c1 <- create_graph(inner_A_c1, vertex_colors)

graph_A_c2 <- create_graph(inner_A_c2, vertex_colors)

graph_A_c3 <- create_graph(inner_A_c3, vertex_colors)

png("inner_graph_A_c1.png",width=3000,height=3000,res=300)
plot(graph_A_c1[[1]])
dev.off()

png("inner_graph_A_c2.png",width=3000,height=3000,res=300)
plot(graph_A_c2[[1]])
dev.off()

png("inner_graph_A_c3.png",width=3100,height=3100,res=300)
plot(graph_A_c3[[1]])
dev.off()

# Subtype top signatures
c1_sigs <- generate_top_signatures(inner_A_c1)

wav_c1_average <- (inner_A_c1_approx + inner_A_c1_h1 + inner_A_c1_h2 + inner_A_c1_h3) / 4
rownames(wav_c1_average) <- rownames(inner_A_c1)
colnames(wav_c1_average) <- colnames(inner_A_c1)
wav_c1_sigs <- generate_top_signatures(wav_c1_average)

c2_sigs <- generate_top_signatures(inner_A_c2)

wav_c2_average <- (inner_A_c2_approx + inner_A_c2_h1 + inner_A_c2_h2 + inner_A_c2_h3) / 4
rownames(wav_c2_average) <- rownames(inner_A_c2)
colnames(wav_c2_average) <- colnames(inner_A_c2)
wav_c2_sigs <- generate_top_signatures(wav_c2_average)


c3_sigs <- generate_top_signatures(inner_A_c3)

wav_c3_average <- (inner_A_c3_approx + inner_A_c3_h1 + inner_A_c3_h2 + inner_A_c3_h3) / 4
rownames(wav_c3_average) <- rownames(inner_A_c3)
colnames(wav_c3_average) <- colnames(inner_A_c3)
wav_c3_sigs <- generate_top_signatures(wav_c3_average)

write.csv(A_sigs[[1]],"A_top_pos_interactions.csv", row.names = FALSE)
write.csv(wav_A_sigs[[1]],"A_top_wav_interactions.csv", row.names = FALSE)
write.csv(S_sigs[[1]],"S_top_pos_interactions.csv", row.names = FALSE)
write.csv(wav_S_sigs[[1]],"S_top_wav_interactions.csv", row.names = FALSE)
write.csv(c1_sigs[[1]],"c1_top_pos_interactions.csv", row.names = FALSE)
write.csv(wav_c1_sigs[[1]],"c1_top_wav_interactions.csv", row.names = FALSE)
write.csv(c2_sigs[[1]],"c2_top_pos_interactions.csv", row.names = FALSE)
write.csv(wav_c2_sigs[[1]],"c2_top_wav_interactions.csv", row.names = FALSE)
write.csv(c3_sigs[[1]],"c3_top_pos_interactions.csv", row.names = FALSE)
write.csv(wav_c3_sigs[[1]],"c3_top_wav_interactions.csv", row.names = FALSE)


