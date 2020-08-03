### MILO FUNCTIONS ###
# Just doing a provisional collection of functions to source in notebooks, will refactor once the class and package is ready
library(edgeR)

countCells <- function(graph, meta, vertex.list, random.vertices, sample.column='Sample'){
  count.matrix <- matrix(0L, ncol=length(unique(meta[, sample.column, drop=TRUE])), nrow=length(vertex.list))
  colnames(count.matrix) <- unique(meta[, sample.column, drop=TRUE])
  
  for(x in seq_along(1:length(vertex.list))){
    v.x <- vertex.list[[x]]
    for(i in seq_along(1:length(unique(meta[, sample.column, drop=TRUE])))){
      i.s <- unique(meta[, sample.column, drop=TRUE])[i]
      i.s.vertices <- intersect(v.x, meta[meta[, sample.column, drop=TRUE] == i.s, ]$Vertex)
      
      count.matrix[x, i] <- length(i.s.vertices)
    }
  }
  rownames(count.matrix) <- rownames(meta)[as.numeric(random.vertices)]
  return(count.matrix)
}

### Testing ###
quant_neighbourhood <- function(graph, meta, sample.column='Sample', sample.vertices=0.25){
  # Count conditions
  random.vertices <- sample(V(graph), size=floor(sample.vertices*length(V(graph))))
  vertex.list <- sapply(1:length(random.vertices), FUN=function(X) neighbors(graph, v=random.vertices[X]))
  
  
  count.matrix <- matrix(0L, ncol=length(unique(meta[, sample.column, drop=TRUE])), nrow=length(vertex.list))
  colnames(count.matrix) <- unique(meta[, sample.column, drop=TRUE])
  
  for(x in seq_along(1:length(vertex.list))){
    v.x <- vertex.list[[x]]
    for(i in seq_along(1:length(unique(meta[, sample.column, drop=TRUE])))){
      i.s <- unique(meta[, sample.column, drop=TRUE])[i]
      i.s.vertices <- intersect(v.x, meta[meta[, sample.column, drop=TRUE] == i.s, ]$Vertex)
      
      count.matrix[x, i] <- length(i.s.vertices)
    }
  }
  rownames(count.matrix) <- random.vertices
  return(list(count.matrix, vertex.list, random.vertices))
}

testQLF <- function(graph, sim2.counts, sim2.model, connectivity='edge', pca=NULL){
  sim2.dge <- DGEList(sim2.counts[, rownames(sim2.model)], lib.size=log(colSums(sim2.counts)))
  sim2.dge <- estimateDisp(sim2.dge, sim2.model)
  sim2.fit <- glmQLFit(sim2.dge, sim2.model, robust=TRUE)
  sim2.contrast <- makeContrasts(ConditionA - ConditionB, levels=sim2.model)
  sim2.res <- glmQLFTest(sim2.fit, contrast=sim2.contrast)
  sim2.res <- as.data.frame(topTags(glmQLFTest(sim2.fit, coef=1), sort.by='none', n=Inf))
  sim2.res$Sig <- as.factor(as.numeric(sim2.res$PValue <= 0.05))
  sim2.res$Neighbourhood <- as.numeric(rownames(sim2.res))
  
  sim2.spatialfdr <- graph_spatialFDR(neighborhoods=vertex.list, graph=graph, connectivity=connectivity, pvalues=sim2.res$PValue, pca = pca)
  return(list(res=sim2.res, spFDR=sim2.spatialfdr))
}


graph_spatialFDR <- function(neighborhoods, graph, pvalues, connectivity='vertex', pca=NULL){
  # input a set of neighborhoods as a list of graph vertices
  # the input graph and the unadjusted GLM p-values
  #' neighborhoods: list of vertices and their respective neighborhoods
  #' graph: input kNN graph
  #' pvalues: a vector of pvalues in the same order as the neighborhood indices
  #' connectivity: character - edge or vertex to calculate neighborhood connectivity or distance to use average Euclidean distance
  #' pca: matrix of PCs to calculate Euclidean distances, only required when connectivity == distance
  # Discarding NA pvalues.
  haspval <- !is.na(pvalues)
  if (!all(haspval)) {
    coords <- coords[haspval, , drop=FALSE]
    pvalues <- pvalues[haspval]
  }
  # define the subgraph for each neighborhood then calculate the vertex connectivity for each
  # this latter computation is quite slow - can it be sped up?
  subgraphs <- lapply(1:length(neighborhoods[haspval]),
                      FUN=function(X) induced_subgraph(graph, neighborhoods[haspval][[X]]))
  # now loop over these sub-graphs to calculate the connectivity - this seems a little slow...
  if(connectivity == "vertex"){
    connect <- lapply(subgraphs, FUN=function(EG) vertex_connectivity(EG))
  } else if(connectivity == "edge"){
    connect <- lapply(subgraphs, FUN=function(EG) edge_connectivity(EG))
  } else if(connectivity == "distance"){
    if(!is.null(pca)){
      connect <- lapply(1:length(neighborhoods[haspval]),
                        FUN=function(PG) {
                          x.pcs <- pca[neighborhoods[haspval][[PG]], ]
                          x.euclid <- as.matrix(dist(x.pcs))
                          x.distdens <- 1/mean(x.euclid[lower.tri(x.euclid, diag=FALSE)])
                          return(x.distdens)})
    } else{
      errorCondition("A matrix of PCs is required to calculate distances")  
    }
  }else{
    errorCondition("connectivity option not recognised - must be either edge or vertex")
  }
  # use 1/connectivity as the weighting for the weighted BH adjustment from Cydar
  w <- 1/unlist(connect)
  w[is.infinite(w)] <- 0
  # Computing a density-weighted q-value.
  o <- order(pvalues)
  pvalues <- pvalues[o]
  w <- w[o]
  adjp <- numeric(length(o))
  adjp[o] <- rev(cummin(rev(sum(w)*pvalues/cumsum(w))))
  adjp <- pmin(adjp, 1)
  if (!all(haspval)) {
    refp <- rep(NA_real_, length(haspval))
    refp[haspval] <- adjp
    adjp <- refp
  }
  return(adjp)
}




### Neighborhood sampling functions ###

sample_nh <- function(milo.sce, 
         graph, ## This will be a slot in the milo class 
         k_param,
         sampling_mode="refined"
         ){
  random_vertices <- 
  if (sampling.mode=="random") {
    sampled_vertices <- random_vertices
  } else if (sampling.mode=="refined") {
    vertex.knn <- BiocNeighbors::findKNN(X=X_pca, k=k_param, subset=as.vector(random_vertices))
    sampled_vertices <- V(graph)[sapply(1:nrow(vertex.knn$index), function(i) refine_vertex(vertex.knn, i, X_pca))]
  }
  sampled_vertices <- unique(sampled_vertices)
  nh_list <- sapply(1:length(sampled_vertices), FUN=function(X) neighbors(graph, v=sampled_vertices[X]))
  nh_list <- setNames(nh_list, sampled_vertices)
  return(nh_list)  
}

refine_vertex <- function(vertex.knn, v.ix, X_pca){
  # vertex.knn: KNN graph for randomly sampled points (output of BiocNeighbors::findKNN)
  # v.ix: index of vertex to refine in vertex.knn
  ## Calculate median profile of KNNs of vertex
  v.med <- apply(X_pca[vertex.knn$index[v.ix,],], 2, median)
  ## Find the closest point to the median and sample
  refined.vertex <- BiocNeighbors::findKNN(rbind(v.med, X_pca), subset=1, k=1)[["index"]][1] - 1 ## -1 to remove the median
  return(refined.vertex)
}



