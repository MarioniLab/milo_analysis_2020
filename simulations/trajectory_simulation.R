#! /usr/bin/env Rscript

library(ggplot2)
library(igraph)
library(ggthemes)
library(ggsci)
library(umap)
library(reshape2)
library(SingleCellExperiment)
library(scran)
library(scater)
library(irlba)
library(mvtnorm)
library(Rfast)
library(miloR)
library(dyntoy)
library(optparse)
library(dplyr)

parser <- OptionParser()

parser <- add_option(parser, c("-n", "--ncells"), type="numeric",
                     help="The total number of cells to simulate")

parser <- add_option(parser, c("-c", "--milestones"), type="numeric",
                     help="The number of milestones to simulate")

parser <- add_option(parser, c("-p", "--props"), type="character",
                     help="A comma-separated list of proportions, one per milestone")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Output directory")

parser <- add_option(parser, c("-x", "--suffix"), type="character",
                     help="Suffix to add to filename to make unique between simulation runs")

opt <- parse_args(parser)

# this is the function to simulate data
simulate_linear_trajectory <- function(n.milestones=3, total.size=100){
    dataset <- generate_dataset(
        model = model_linear(num_milestones = n.milestones),
        num_cells = total.size,
        num_features = 2000
        )
    
    ## Build SingleCellExperiment object
    cnts <- t(dataset$counts)
    coldata <- data.frame(row.names = colnames(cnts), dataset$prior_information$groups_id)
     
    ## Dimensionality reduction
    pca <- prcomp_irlba(dataset$expression, n=50, scale.=TRUE, center=TRUE)
    sce <- SingleCellExperiment(assays=list(counts=cnts, logcounts=t(dataset$expression)), colData=coldata,
                                reducedDims=list("PCA"=pca$x))
    mylo <- Milo(sce)

    # add condition information
    branches <- dataset$prior_information$groups_id
    print(head(branches))
    print(head(rownames(mylo)))
    coldata_df <- data.frame(cell_id = colnames(mylo))
    coldata_df <- left_join(coldata_df, branches)
    print(head(coldata_df))
    ## Simulate biological condition
    n_groups <- length(unique(branches$group_id))
    p_vec <- as.numeric(unlist(strsplit(opt$props, split=",", fixed=TRUE)))
    if(length(p_vec) == 1){
        p_vec <- rep(p_vec, opt$milestones)
    } else if(length(p_vec) != opt$milestones){
        stop("Condition proportions must == length(milestones)")
    }
    print(p_vec)
    a.cells <- c()
    for (i in 1:n_groups) {
      g <- paste0("M",i)
      p <- p_vec[i]

      print(g)
      print(sum(coldata_df$group_id==g))
      print(floor(sum(coldata_df$group_id==g)*p))
      print(length(coldata_df$cell_id[coldata_df$group_id==g]))

      m.A <- sample(coldata_df$cell_id[coldata_df$group_id==g], 
                size=floor(sum(coldata_df$group_id==g)*p))
      a.cells <- c(a.cells, m.A)
   }
   
   coldata_df <- coldata_df %>% dplyr::mutate(Condition = ifelse(cell_id %in% a.cells, 'A', 'B'))

   ## Simulate replicates
   coldata_df <- coldata_df %>%
     group_by(group_id) %>%
     dplyr::mutate(Replicate=c(rep("R1", floor(n()*0.3)), 
                               rep("R2", floor(n()*0.3)), 
                               rep("R3", n() - 2*(floor(n()*0.3))))
  )
  coldata_df$Sample <- paste(coldata_df$Condition, coldata_df$Replicate, sep="_")

  return(list("mylo"=mylo, "meta"=as.data.frame(coldata_df)))
}

message(paste0("Simulating ", opt$ncells, " cells in a single trajectory across ", opt$milestones, " milestones"))
x.traj <- simulate_linear_trajectory(n.milestones=opt$milestones, total.size=opt$ncells)

ofile <- paste(opt$output, paste("Trajectory", paste0("Ncells", opt$ncells), opt$suffix, "simpleSim.RDS", sep="_"), sep="/")
message(paste0("Saving simulated data to: ", ofile))
saveRDS(x.traj, file=ofile)