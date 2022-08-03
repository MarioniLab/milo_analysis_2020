#! /usr/bin/env Rscript

library(dyngen)
library(dplyr)
library(rlang)
library(SingleCellExperiment)
library(Matrix)
library(umap)
library(scran)
library(irlba)
library(scater)

backbone <- backbone_branching(
  num_modifications = 2,
  min_degree = 3,
  max_degree = 3
)

init <- initialise_model(
  backbone = backbone,
  num_cells = 5500,
  num_tfs = 600,
  num_targets = 3650,
  num_hks = 350,
  simulation_params = simulation_default(census_interval = 4, ssa_algorithm = ssa_etl(tau = 300 / 3600)),
  download_cache_dir = "~/.cache/dyngen/",
  num_cores = 6
)
branch.sim <- generate_dataset(init, make_plots = TRUE)

saveRDS(branch.sim, "BranchingTrajectory.RDS")

genes <- colnames(branch.sim$dataset$counts)
names(genes) <- NULL

branch.counts <- as(t(branch.sim$dataset$counts), "dgCMatrix")
branch.log <- as(t(branch.sim$dataset$expression), "dgCMatrix")
rownames(branch.counts) <- genes
rownames(branch.log) <- genes

branch.sce <- SingleCellExperiment(assay=list(counts=branch.counts, logcounts=branch.log))
reducedDim(branch.sce, "Manifold") <- branch.sim$dataset$dimred

# assign milestones to cells in different regions
# the starting milestones are in the progressions slot column 'from'
branch.milestones <- as.data.frame(branch.sim$dataset$progressions[, c("cell_id", "from")])
rownames(branch.milestones) <- branch.milestones$cell_id
branch.milestones <- branch.milestones[colnames(branch.sce), ]
colData(branch.sce)$Milestone <- branch.milestones$from

hvgs <- modelGeneVar(branch.sce)
hvgs$FDR[is.na(hvgs$FDR)] <- 1
# top.hvgs <- rownames(hvgs[order(hvgs$FDR, decreasing=FALSE), ])[1:100]
top.hvgs <- hvgs$FDR < 0.05

branch.pca <- prcomp_irlba(t(logcounts(branch.sce[top.hvgs, ])), n=50)
reducedDim(branch.sce, "PCA") <- branch.pca$x

set.seed(42)
branch.umap <- umap(branch.pca$x[, c(1:30)],
                    n_neighbours=21,
                    init='random', min_dist=0.3)$layout
reducedDim(branch.sce, "UMAP") <- branch.umap

saveRDS(branch.sce, file="BranchingTrajectory_SCE.RDS")








