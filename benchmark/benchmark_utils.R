### BENCHMARKING FUNCTIONS ###

suppressPackageStartupMessages({
library(SingleCellExperiment)
library(DAseq)
library(miloR)
library(tibble)
library(dplyr)
library(igraph)
library(cydar)
library(pdist)
  library(reshape2)
  })

# # ## Set-up reticulate 4 MELD
# reticulate::use_condaenv("emma_env", required=TRUE)
# library(reticulate) ## development version of reticulate, or numba use breaks C stack

### SYNTHETIC LABELS ###

.find_centroid <- function(X_emb, cluster_membership){
  cl.ixs <- split(1:nrow(X_emb), cluster_membership)  
  centroid_emb <- sapply(cl.ixs, function(x) colMeans(X_emb[x, , drop=FALSE]))
  centroid_emb
}

.member_weight <- function(x, centroid_dist, m=2){
  # centroid_dist <- pdist(t(x), t(centroid_emb))@dist
  w_memb <- sapply(centroid_dist, function(x) 1/sum(x/centroid_dist)^(2/(m-1)))
}

.scale_to_range <- function(x, min=1, max=10){
  ((x - min(x))/(max(x)-min(x)))*(max-min) + min
}


.logit <- function(x, a=1){
  1/(1+exp(- a * x))
}

# Creates random differentially expressed regions over a dataset for benchmarking.
add_synthetic_labels_pop <- function(sce, # SingleCellExperiment obj
                                     pop, pop_column="celltype",
                                     pop_enr = 0.7,
                                     redDim='pca.corrected', # embedding to use to simulate differential abundance
                                     n_conditions=2, # number of conditions to simulate
                                     n_replicates=3, # number of replicates per condition
                                     n_batches = 2, # number of technical batches per condition (at least 2 replicates per batch)
                                     condition_balance = 1, # the distribution of cells across conditions (1 = equal)
                                     m=2, # Fuzziness parameter (higher m, more fuzziness)
                                     a_logit=0.5, # logit parameter
                                     cap_enr=NULL,
                                     seed=42){
  
  # pop_sce = sce[,sce[[pop_column]]==pop]
  
  set.seed(seed)
  conditions = paste0("Condition", 1:n_conditions)
  
  X_emb = reducedDim(sce, redDim)
  
  ## Find cluster center
  cluster_membership = sce[[pop_column]]
  centroid_emb <- .find_centroid(X_emb, cluster_membership)

  ## Assign weight to each cell for each cluster center
  centroid_dist <- pdist(X_emb, t(centroid_emb))
  centroid_dist <- as.matrix(centroid_dist)
  
  w <- sapply(1:ncol(centroid_dist),  function(j) 
    sapply(1:nrow(centroid_dist), function(i) 
      1/sum(centroid_dist[i,j]/centroid_dist[i,])^(2/(m-1))
    ) 
  )
  colnames(w) <- colnames(centroid_emb)
  rownames(w) <- rownames(X_emb)
  w <- apply(scale(w), 2, .logit, a=a_logit)
  ## Normalize weights from enr_score to 0.5
  enr_scores <- rep(0.5, ncol(w)) ## Generate enrichment prob for each cluster
  # enr_scores <- runif(ncol(w)) ## Generate enrichment prob for each cluster
  names(enr_scores) <- colnames(w)
  if(length(pop_enr) == length(pop)){
    enr_scores[pop] <- pop_enr
  } else{
    # assume all pops have the same enrichment
    pop_enr <- rep(pop_enr, length(pop))
    enr_scores[pop] <- pop_enr
  }
  
  
  # altering the baseline probability can induce a skew towards a condition across _all_ cells
  enr_prob <- sapply(1:ncol(w), function(i) .scale_to_range(w[,i], min=0.5*condition_balance,
                                                            max=enr_scores[i]))
  colnames(enr_prob) <- colnames(centroid_emb)
  
  # need to integrate over these to get the condition probabilities
  # need to set relevant pops only, force the others to ~0.5
  prob_matrix <- enr_prob[,pop]
  if(is(prob_matrix, "matrix")){
    cond_probability <- rowMeans(prob_matrix)
    for(x in seq_along(pop)){
      cond_probability[sce[[pop_column]] == pop[x]] <- prob_matrix[sce[[pop_column]] == pop[x], pop[x]]
    }
  } else{
    cond_probability <- prob_matrix
  }
  
  ## Cap probabilities (to have the same number of DA cells w different maximum Fold Change)
  if (!is.null(cap_enr)) {
    cond_probability <- ifelse(cond_probability > cap_enr, cap_enr, cond_probability)
  }
  sim3$Condition1_prob <- ifelse(sim3$Condition1_prob > 0.8, 0.8, sim3$Condition1_prob)
  
  cond_probability = cbind(cond_probability, 1 - cond_probability)
  colnames(cond_probability) = conditions
  
  # Generate labels for condition and replicates
  synth_labels <- sapply(1:nrow(cond_probability),  function(i) sample(colnames(cond_probability), size = 1, prob = cond_probability[i,]))
  
  replicates <- paste0("R", 1:n_replicates)
  batches <- sample(paste0("B", rep(1:n_batches, each=n_replicates)))
  
  synth_samples <- paste0(synth_labels, "_", replicates)
  if(n_batches > 1){
   names(batches) <- sort(unique(synth_samples))
  } else{
    batches <- rep("B1", length(unique(synth_samples)))
    names(batches) <- unique(synth_samples)
  }
  synth_batches <- batches[synth_samples]
  
  # Add synthetic labels and probabilities to colData
  colData(sce)[["synth_labels"]] <- synth_labels
  # colData(sce)[["synth_replicates"]] <- synth_replicates
  colData(sce)[["synth_samples"]] <- synth_samples
  colData(sce)[["synth_batches"]] <- synth_batches
  colnames(cond_probability) <- paste0(colnames(cond_probability), "_prob")
  colData(sce)[colnames(cond_probability)] <- cond_probability
  return(sce)
}

## Group true DA cells in clusters (to test coverage of DA methods)
cluster_synthetic_labels <- function(embryo_sce, graph, min_cl_size=5){
  adj <- get.adjacency(graph)
  
  ## Retain DA cells
  da.adj <- adj[embryo_sce$true_labels!='NotDA',embryo_sce$true_labels!='NotDA']
  
  ## REmove edges between cells with discodant LFC sign
  da.cells.mat <- sapply(unique(embryo_sce$true_labels), function(x) as.numeric(embryo_sce$true_labels==x))
  da.cells.mat <- da.cells.mat[embryo_sce$true_labels!='NotDA',c("NegLFC", "PosLFC")]
  concord.sign.adj <- tcrossprod(da.cells.mat[,c("NegLFC", "PosLFC")], da.cells.mat[,c("NegLFC", "PosLFC")])
  concord.sign.adj <- as(concord.sign.adj, 'sparseMatrix')
  da.adj[concord.sign.adj==0] <- 0
  
  ## Cluster DA cells
  da.graph <- graph_from_adjacency_matrix(da.adj, mode = 'undirected')
  clust <- igraph::cluster_louvain(da.graph)
  embryo_sce$true_DA_clust <- rep(NA, length(embryo_sce$true_labels))
  embryo_sce$true_DA_clust[embryo_sce$true_labels != "NotDA"] <- clust$membership
  
  ## Remove singletons (or less than min_cl_size cells)
  embryo_sce$true_DA_clust[embryo_sce$true_DA_clust %in% which(table(clust$membership) < min_cl_size)] <- NA
  
  embryo_sce
}

### SYNTHETIC BATCH EFFECT ###

add_batch_effect <- function(embryo_sce, batch_col="synth_samples", norm_sd=0.5){
  cellids_sample <- split(embryo_sce$cell, embryo_sce[[batch_col]])
  X_pca <- reducedDim(embryo_sce, "pca.corrected")
  X_pca_batch <- X_pca

  for (b in names(cellids_sample)){
    batch_effect <- rnorm(ncol(X_pca), mean=0, sd = norm_sd)
    X_pca_batch[cellids_sample[[b]],] <- t(apply(X_pca_batch[cellids_sample[[b]],], 1, function(x) x + batch_effect))
  }
  
  reducedDim(embryo_sce, "pca_batch") <- X_pca_batch
  embryo_sce  
}

### METHODS ###

## Milo

run_milo <- function(sce, condition_col, sample_col, reduced.dim="PCA",
                     k=15, d=30, prop=0.1, returnMilo = TRUE,
                     batch_col=NULL){
  ## Make design matrix
  design_df <- as.tibble(colData(sce)[c(sample_col, condition_col, batch_col)]) %>%
    distinct() %>%
    column_to_rownames(sample_col)
  if (is.null(batch_col)) {
    design <- formula(paste('~', condition_col, collapse = ' '))  
  } else {
    design <- formula(paste('~', batch_col, "+", condition_col, collapse = ' '))
  }
  
  ## Build graph neighbourhoods
  milo <- Milo(sce)
  milo <- buildGraph(milo, k=k, d=d, reduced.dim = reduced.dim)
  milo <- makeNhoods(milo, prop = prop, k=k, d=d, reduced_dims = reduced.dim)
  ## Test DA
  milo <- miloR::countCells(milo, meta.data = data.frame(colData(milo)), sample=sample_col)
  milo <- calcNhoodDistance(milo, d=d, reduced.dim = reduced.dim)
  DA_results <- testNhoods(milo, design = design, design.df = design_df)
  if (isTRUE(returnMilo)) {
    return(list(Milo=milo, DAres=DA_results))
  } else {
    DA_results 
  }
}

milo2output <- function(milo, da_res, out_type="continuous", alpha=0.1){
  if (out_type=="continuous") { 
    da.cell.mat <- milo@nhoods %*% da_res$logFC
    da.cell <- da.cell.mat[,1]
  } else {
    da.nhoods <- ifelse(da_res$SpatialFDR < alpha, ifelse(da_res$logFC > 0, "PosLFC", 'NegLFC'), "NotDA")
    da.nhoods.mat <- sapply(unique(da.nhoods), function(x) as.numeric(da.nhoods==x))
    da.cell.mat <- milo@nhoods %*% da.nhoods.mat
    da.cell <- apply(da.cell.mat, 1, function(x) colnames(da.cell.mat)[which.max(x)])
  }
  da.cell
}


## DAseq

run_daseq <- function(sce, k.vec, condition_col, cell_col="cell", reduced.dim = "PCA", d=30){
  condition_vec <- colData(sce)[[condition_col]]
  conds <- unique(condition_vec)
  cell.labels <- sce[[cell_col]]
  
  daseq_res <- getDAcells(X=reducedDim(sce, reduced.dim)[,1:d],
                          cell.labels=cell.labels,
                          labels.1=cell.labels[condition_vec == conds[2]],
                          labels.2=cell.labels[condition_vec == conds[1]],
                          k.vector=k.vec,
                          size=1)
  return(daseq_res)
  }

daseq2output <- function(sce, daseq_res, cell_col="cell", out_type="continuous"){
  if (out_type=="continuous") {
    da.cell <- daseq_res$da.pred
  } else {
    cell.labels <-  sce[[cell_col]]
    da.cell <- rep("NotDA", length(cell.labels))
    da.cell[daseq_res$da.up] <- "PosLFC"
    da.cell[daseq_res$da.down] <- "NegLFC"
  }
  da.cell
}

## MELD

run_meld_reticulate <- function(sce, condition_col, sample_col, reduced.dim="PCA", d=30, k=15) {
  ## Extract input from SingleCellExperiment
  X_red_dim = reducedDim(sce, reduced.dim)[,1:d]
  condition_vec <- colData(sce)[[condition_col]]
  sample_labels <- colData(sce)[[sample_col]]
  conditions <- sort(unique(condition_vec))
  ## Run MELD
  reticulate::source_python("./run_meld.py")
  py$run_meld(X_red_dim, sample_labels, conditions, k=k)
}

meld2output <- function(meld_res, likelihood_cutoff = 0.6, out_type="continuous"){
  conds <- colnames(meld_res)
  if (out_type=="continuous") {
    da.cell <- meld_res[conds[1]][,1]
  } else {
  da.cell <- rep("NotDA", nrow(meld_res))
  da.cell[meld_res[conds[1]] > likelihood_cutoff] <- "PosLFC"
  da.cell[meld_res[conds[2]] > likelihood_cutoff] <- "NegLFC"
  }
  da.cell
}

## Cydar


run_cydar <- function(sce, condition_col="synth_labels",
                      sample_col="synth_samples",
                      reduced.dim="pca.corrected",
                      d=20,
                      batch_col=NULL,
                      alpha=0.1,
                      tol=1.0,
                      downsample=10,
                      returnCd=TRUE){
  ## Make design matrix
  design_df <- as.tibble(colData(sce)[c(sample_col, condition_col, batch_col)]) %>%
    distinct() %>%
    column_to_rownames(sample_col)
  if (is.null(batch_col)) {
    design <- formula(paste('~', condition_col, collapse = ' '))  
  } else {
    design <- formula(paste('~', batch_col, "+", condition_col, collapse = ' '))
  }
  
  ## Make list for each sample
  sample_ls <- split(1:ncol(sce), sce[[sample_col]])
  processed.exprs <- lapply(sample_ls, function(s) reducedDim(sce[,s], reduced.dim)[,1:d])
  cd <- prepareCellData(processed.exprs)
  ## Count cells in hyperspheres
  cd <- cydar::countCells(cd, tol=tol, filter=1, downsample=downsample)
  # do DA testing with edgeR
  cd.dge <- DGEList(assay(cd), lib.size=cd$totals)
  
  sim.design <- model.matrix(design, data=design_df)[colnames(cd.dge),]
  sim.dge <- estimateDisp(cd.dge, sim.design)
  sim.fit <- glmQLFit(sim.dge, sim.design)
  sim.res <- glmQLFTest(sim.fit, coef=2)
  
  # control the spatial FDR
  cydar.res <- sim.res$table
  cydar.res$SpatialFDR <- spatialFDR(intensities(cd), sim.res$table$PValue)
  is.sig <- cydar.res$SpatialFDR <= alpha
  if (returnCd) {
    list(Cd=cd, DAres=cydar.res)
  } else {
    cydar.res
  }
}

cydar2output <- function(sce, cd, da_res, out_type="continuous", alpha=0.1, sample_col="synth_samples", reduced.dim="pca.corrected", d=30){
  nhs <- lapply(cellAssignments(cd), function(hs) as.vector(hs))
  # hs_mat <- sapply(nhs, function(nh) ifelse(1:sum(cd@colData$totals) %in% nh, 1, 0))
  ## Recover cell ids
  ordered.cells <- colnames(cellIntensities(cd))
  hs_mat <- sapply(nhs, function(nh) ifelse(seq_along(ordered.cells) %in% nh, 1, 0))
  rownames(hs_mat) <- ordered.cells
  colnames(hs_mat) <- rownames(da_res)
  if (out_type=="continuous") { 
    da.cell.mat <- hs_mat %*% da_res$logFC
    da.cell <- da.cell.mat[,1]
  } else {
    da.hs <- ifelse(da_res$SpatialFDR < alpha, ifelse(da_res$logFC > 0, "PosLFC", 'NegLFC'), "NotDA")
    da.hs.mat <- sapply(unique(da.hs), function(x) as.numeric(da.hs==x))
    da.cell.mat <- hs_mat %*% da.hs.mat
    da.cell <- apply(da.cell.mat, 1, function(x) colnames(da.cell.mat)[which.max(x)])
  }
  da.cell[colnames(sce)]
}

## Louvain clustering

run_louvain <- function(sce, condition_col, sample_col, k=15, d=30, reduced.dim="PCA", batch_col=NULL){
  ## Make design matrix
  design_df <- as.tibble(colData(sce)[c(sample_col, condition_col, batch_col)]) %>%
    distinct() %>%
    dplyr::rename(sample=sample_col)
  if (is.null(batch_col)) {
    design <- formula(paste('Freq ~', condition_col, "+ offset(log(N_s))", collapse = ' '))  
  } else {
    design <- formula(paste('Freq ~', batch_col, "+", condition_col, "+ offset(log(N_s))", collapse = ' '))
  }
  ## Louvain clustering
  X_red_dim = reducedDim(sce, reduced.dim)[,1:d]
  sce.graph <- buildKNNGraph(t(X_red_dim), k=k)
  louvain.clust <- cluster_louvain(sce.graph)
  louvain.clust.ids <- membership(louvain.clust)
  
  condition_vec <- colData(sce)[[condition_col]]
  sample_labels <- colData(sce)[[sample_col]]
  clust.df <- data.frame("cell_id"=colnames(sce), "Louvain.Clust"=as.character(louvain.clust.ids))
  clust.df$Sample <- sample_labels
  clust.df$Condition <- condition_vec
  
  louvain.count <- table(clust.df$Louvain.Clust, clust.df$Sample)
  attributes(louvain.count)$class <- "matrix"
  
  df <- melt(louvain.count, varnames=c("cluster", "sample"),  value.name="Freq") %>%
    mutate(cluster=factor(cluster)) %>%
    left_join(design_df, by="sample") %>%
    group_by(sample) %>%
    mutate(N_s=sum(Freq)) %>%
    ungroup() %>%
    group_by(cluster) %>%
    do(model=glm(design, data=., family="poisson")) 
  
  res_df <- t(sapply(df$model, function(x) summary(x)$coefficients[nrow(summary(x)$coefficients),]))
  colnames(res_df) <- c("logFC","Std. Error", "z value",    "Pval" )
  louvain.res <- cbind(df, res_df) %>%
    mutate(FDR=p.adjust(Pval, method = "BH"))
  rownames(louvain.res) <- louvain.res$cluster
  
  clust.df$logFC <- louvain.res[clust.df$Louvain.Clust, 'logFC']
  clust.df$FDR <- louvain.res[clust.df$Louvain.Clust, 'FDR']
  return(clust.df)
  }

louvain2output <- function(louvain_res, out_type="continuous", alpha=0.1, lfc_threshold=0){
  if (out_type=="continuous") {
    da.cell <- louvain_res$logFC
  } else {
    da.cell <- ifelse(louvain_res$FDR < alpha & abs(louvain_res$logFC) > lfc_threshold, ifelse(louvain_res$logFC > 0, "PosLFC", 'NegLFC'), "NotDA")
  }
  da.cell
}

### RUN BENCHMARK ON SYNTHETIC LABELS ###

runDA <- function(sce, coldata, X_pca, 
                  method, 
                  condition_col='synth_labels', 
                  sample_col="synth_samples",
                  params = list(milo = list(k=15),
                                meld = list(k=15),
                                daseq = list(k.vec=c(10,20,30, 40)),
                                louvain = list(k=15),
                                cydar = list(tol=5, downsample=3)
                  ), d=30, out_type="label"
                  ){
  ## Check that method name is in params
  if (!method %in% names(params)) {
    stop(paste("Specify parameters for method", method))
  }
  ## Check valid method
  if (!method %in% c("milo", "milo_batch", "louvain", 'daseq', "louvain_batch", 'cydar', 'cydar_batch')) {
    stop("Unrecognized method")
  }
  
  ## Add reduced dim + coldata to sce
  colData(sce) <- DataFrame(coldata)
  reducedDim(sce, "pca_batch") <- as.matrix(X_pca)
  
  ## Run method
  ## Run milo
  if (method=="milo") {
    milo_res <- run_milo(sce, condition_col=condition_col, sample_col=sample_col,
                         reduced.dim = "pca_batch", d=d, k=params$milo$k)
    out <- milo2output(milo_res$Milo, milo_res$DAres, out_type = out_type)
  } else if (method == "milo_batch"){
    ## Run milo controlling for batch
    milo_batch_res <- run_milo(sce, condition_col=condition_col, sample_col=sample_col,
                               reduced.dim = "pca_batch", d=d, k=params$milo$k, batch_col = "synth_batches")
    out <- milo2output(milo_batch_res$Milo, milo_batch_res$DAres, out_type = out_type)
  } else if (method == "cydar"){
    cydar_res <- run_cydar(sce, condition_col=condition_col, sample_col=sample_col,
                               reduced.dim = "pca_batch", d=d, tol=params$cydar$tol,
                           downsample=params$cydar$downsample)
    out <- cydar2output(sce, cydar_res$Cd, cydar_res$DAres, out_type = out_type)
  } else if (method == "cydar_batch"){
    cydar_batch_res <- run_cydar(sce, condition_col=condition_col, sample_col=sample_col,
                           reduced.dim = "pca_batch", d=d, tol=params$cydar$tol,
                           downsample=params$cydar$downsample, batch_col = "synth_batches")
    out <- cydar2output(cydar_batch_res$Cd, cydar_batch_res$DAres, out_type = out_type)
  } else if (method == "daseq"){
    ## Run DAseq
    daseq_res <- run_daseq(sce, k.vec=params$daseq$k.vec, condition_col, 
                           reduced.dim = "pca_batch", d=d)
    out <- daseq2output(sce, daseq_res, out_type = out_type)
  } else if (method == "louvain"){
    ## Run louvain
    louvain_res <- run_louvain(sce, condition_col=condition_col, sample_col=sample_col,
                               reduced.dim = "pca_batch", d=d, k=params$louvain$k)
    out <- louvain2output(louvain_res, out_type = out_type)
  } else if (method == "louvain_batch"){
    ## Run louvain
    louvain_batch_res <- run_louvain(sce, condition_col=condition_col, sample_col=sample_col,
                               reduced.dim = "pca_batch", d=d, k=params$louvain$k, batch_col="synth_batches")
    out <- louvain2output(louvain_batch_res, out_type = out_type)
  }
  
  ## Save results
  bm <- data.frame(bm_out=out)
  if (length(out) < ncol(sce)) {
    return(out)
  }
  bm$true_prob <- sce$Condition2_prob 
  bm$true <- sce$true_labels
  if (!is.null(sce$true_DA_clust)) {
    bm$true_clust <- sce$true_DA_clust
  }
  long_bm <- pivot_longer(bm, cols = bm_out, names_to='method', values_to="pred")
  long_bm[["method"]] <- method
  return(long_bm)
}


calculate_outcome <- function(long_bm){
  long_bm <- long_bm %>%
    mutate(outcome=case_when(true==pred & pred!="NotDA" ~ 'TP',
                             true!=pred & pred!="NotDA" ~ 'FP',
                             true!=pred & pred=="NotDA" ~ 'FN',
                             true==pred & pred=="NotDA"  ~ "TN"
    )) %>%
    group_by(method, outcome) %>%
    summarise(n=n()) %>%
    pivot_wider(id_cols=method, names_from=outcome, values_from=n, values_fill=0) 
  
  check_cols <- c("TP","FP","FN","TN") %in% colnames(long_bm)
  if (any(!check_cols)) {
    add_cols <- c("TP","FP","FN","TN")[!check_cols]
    for (col in add_cols) {
      long_bm[[col]] <- rep(0, nrow(long_bm))
      }
  }
  
  long_bm %>%
    mutate(TPR=TP/(TP+FN), FPR=FP/(FP+TN), TNR=TN/(TN+FP), FNR = FN/(FN+TP),
           FDR = FP/(TP+FP),
           Precision = TP/(TP+FP),
           Power = 1 - FNR,
           Accuracy = (TP + TN)/(TP + TN + FP + FN)
    )
}

## --- old functions --- ##

benchmark_da <- function(sce, condition_col='synth_labels', 
                         sample_col="synth_samples",
                         red_dim="pca.corrected",
                         params = list(milo = list(k=15),
                                       meld = list(k=15),
                                       daseq = list(k.vec=c(10,20,30, 40)),
                                       louvain = list(k=15),
                                       cydar = list(tol=5, downsample=3)
                         ),
                         d=30, out_type = "continuous"){
  ## Run milo
  milo_res <- run_milo(sce, condition_col=condition_col, sample_col=sample_col,
                       reduced.dim = red_dim, d=d, k=params$milo$k)
  milo_out <- milo2output(milo_res$Milo, milo_res$DAres, out_type = out_type)
  ## Run milo controlling for batch
  milo_batch_res <- run_milo(sce, condition_col=condition_col, sample_col=sample_col,
                             reduced.dim = red_dim, d=d, k=params$milo$k, batch_col = "synth_batches")
  milo_batch_out <- milo2output(milo_batch_res$Milo, milo_batch_res$DAres, out_type = out_type)
  
  ## Run DAseq
  daseq_res <- run_daseq(sce, k.vec=params$daseq$k.vec, condition_col, 
                         reduced.dim = red_dim, d=d)
  daseq_out <- daseq2output(sce, daseq_res, out_type = out_type)
  ## Run MELD
  meld_res <- run_meld_reticulate(sce, condition_col=condition_col, sample_col=sample_col,
                                  reduced.dim = red_dim, d=d, k=params$meld$k)
  meld_out <- meld2output(meld_res, out_type = out_type)
  ## Run louvain
  louvain_res <- run_louvain(sce, condition_col=condition_col, sample_col=sample_col,
                             reduced.dim = red_dim, d=d, k=params$louvain$k)
  louvain_out <- louvain2output(louvain_res, out_type = out_type)
  ## Collect results + true labels
  bm <- data.frame(milo=milo_out, milo_batch=milo_batch_out,
                   daseq=daseq_out, meld=meld_out, louvain=louvain_out)
  bm$true_prob <- sce$Condition2_prob 
  bm$true <- sce$true_labels
  if (!is.null(sce$true_DA_clust)) {
    bm$true_clust <- sce$true_DA_clust
  }
  long_bm <- pivot_longer(bm, cols = c(milo, milo_batch, daseq, meld, louvain), names_to='method', values_to="pred") 
  return(long_bm)
}


# Given a data_embedding, sample a simplex to weight each dimension to get a
# new PDF for each condition
# a: scaling coefficient in logit (lower a --> less extreme probabilities) 
.make_pdf <- function(data_embedding, a=0.2){
  # Create an array of values that sums to 1
  n_components = ncol(data_embedding)
  data_simplex = sort(runif(n = n_components-1))
  data_simplex = c(0, data_simplex, 1)
  data_simplex = diff(data_simplex)
  data_simplex = sample(data_simplex)
  # Weight each embedding component by the simplex weights
  sort_axis = rowSums(data_embedding * data_simplex)
  
  # Pass the weighted components through a logit
  pdf = 1/(1+exp(- a * sort_axis))
  if (sample(c(TRUE, FALSE), 1)){ pdf = 1 - pdf }
  return(pdf)  
}

## Smooth probabilities over KNN graph (to avoid having clusters with opposite sign DA)
.knn_smoothing <- function(graph, cond_probability, k=15, redDim='pca.corrected', d=10){
  # X_red_dim = reducedDim(sce, redDim)[,1:d]
  # graph = buildKNNGraph(t(X_red_dim), k = k)  
  adj = get.adjacency(graph)
  smooth_cond_probability <- (adj %*% cond_probability)/rowSums(adj)
  smooth_cond_probability
}

# Creates random differentially expressed regions over a dataset for benchmarking.
add_synthetic_labels <- function(sce, # SingleCellExperiment obj
                                 knn_graph, # for knn smoothing of probability values
                                 redDim='pca.corrected', # embedding to use to simulate differential abundance
                                 n_conditions=2, # number of conditions to simulate
                                 n_components=10, # number of components of embedding to use
                                 n_replicates=3, # number of replicates per condition
                                 n_batches = 2, # number of technical batches per condition (at least 2 replicates per batch)
                                 seed=42){
  data_embedding = reducedDim(sce, redDim)[,1:n_components]
  set.seed(seed)
  # embedding data must be mean-centered
  data_embedding = t(scale(t(data_embedding), scale=FALSE))
  
  # Randomly flip sign of each embedding dimension
  data_embedding = apply(data_embedding, 2, function(x)  x*sample(c(-1, 1), 1) )
  
  conditions = paste0("Condition", 1:n_conditions)
  cond_probability = sapply(1:(length(conditions)-1), function(x) .make_pdf(data_embedding))
  
  # KNN Smoothing to avoid regions of the graph with opposite labels
  cond_probability = .knn_smoothing(knn_graph, cond_probability, redDim=redDim, d=n_components)
  
  # Normalize to sum to 1 for each cell
  # cond_probability <- t(apply(cond_probability, 1, function(x) x/sum(abs(x))))
  cond_probability = cbind(cond_probability, 1 - rowSums(cond_probability))
  colnames(cond_probability) = conditions
  
  # Generate labels for condition and replicates
  synth_labels <- sapply(1:nrow(cond_probability),  function(i) sample(colnames(cond_probability), size = 1, prob = cond_probability[i,]))
  replicates <- paste0("R", 1:n_replicates)
  batches <- sample(paste0("B", rep(1:n_batches, each=n_replicates)))
  synth_samples <- paste0(synth_labels, "_", replicates)
  names(batches) <- sort(unique(synth_samples))
  synth_batches <- batches[synth_samples]
  
  # Add synthetic labels and probabilities to colData
  colData(sce)[["synth_labels"]] <- synth_labels
  # colData(sce)[["synth_replicates"]] <- synth_replicates
  colData(sce)[["synth_samples"]] <- synth_samples
  colData(sce)[["synth_batches"]] <- synth_batches
  colnames(cond_probability) <- paste0(colnames(cond_probability), "_prob")
  colData(sce)[colnames(cond_probability)] <- cond_probability
  return(sce)
}

# ## Generate pdf along a randomly selected diffusion component
# .diffmap2pdf <- function(sce, ncomponents=5, a=1, b=1, d=30){
#   sce <- runDiffusionMap(sce, dimred="pca.corrected", ncomponents=ncomponents, n_dimred=d)
#   diff_comp = reducedDim(sce, 'DiffusionMap')[,sample(1:ncomponents,1)]
#   ## Convert latent dimension to probability
#   cond_probability = 1/(1+a*exp( - b * scale(diff_comp)))
#   cond_probability
# }
