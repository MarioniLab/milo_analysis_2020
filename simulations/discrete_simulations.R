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
library(optparse)

parser <- OptionParser()

parser <- add_option(parser, c("-n", "--ncells"), type="numeric",
                     help="The total number of cells to simulate")

parser <- add_option(parser, c("-c", "--clusters"), type="numeric",
                     help="The number of clusters to simulate")

parser <- add_option(parser, c("-p", "--props"), type="character",
                     help="A comma-separated list of proportions, one per cluster")

parser <- add_option(parser, c("-l", "--cellsPerCluster"), type="numeric",
       	             help="A comma-separated list of the number of cells per cluster. If not set then assumed to be ~balanced.")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Output directory")

parser <- add_option(parser, c("-x", "--suffix"), type="character",
                     help="Suffix to add to filename to make unique between simulation runs")

opt <- parse_args(parser)


# return a Milo object of the expected size
simulate_discrete_data <- function(n.clusters=3, total.size=100, cells.per.cluster=c(30, 30, 40),
                                   condition.props=c(0.9, 0.1, 0.5)){
    set.seed(42)
    r.n <- 1000
    n.dim <- 50
    if(n.clusters != length(cells.per.cluster)){
        stop("Number of clusters must equal the length of cells.per.cluster")
    }
    
    if(total.size != sum(cells.per.cluster)){
        stop("Total sample size must match the sum of cells.per.cluster")
    }
    
    gex.list <- list()
    for(x in seq_along(1:n.clusters)){
        block.cells <- cells.per.cluster[x]
        # select a set of eigen values for the covariance matrix of each block, say 50 eigenvalues?
        # randomly sample a mean value
        block.mean <- runif(n=1, min=2, max=7)
        block.eigens <- sapply(1:n.dim, FUN=function(X) rexp(n=1, rate=abs(runif(n=1, min=0, max=50))))
        block.eigens <- block.eigens[order(block.eigens)]
        block.p <- qr.Q(qr(matrix(rnorm(block.cells^2, mean=4, sd=0.01), block.cells)))
        block.sigma <- crossprod(block.p*block.eigens, block.p*block.eigens)
        block.gex <- abs(Rfast::rmvnorm(n=r.n, mu=rnorm(n=block.cells, mean=block.mean, sd=0.01), sigma=block.sigma))
        gex.list[[paste0("Block", x)]] <- block.gex
        
    }
    
    sim.gex <- do.call(cbind, gex.list)
    colnames(sim.gex) <- paste0("Cell", 1:ncol(sim.gex))
    rownames(sim.gex) <- paste0("Gene", 1:nrow(sim.gex))
    sim.pca <- prcomp_irlba(t(sim.gex), n=50, scale.=TRUE, center=TRUE)
    
    if(length(condition.props) != length(cells.per.cluster)){
        stop("The length of condition.props must be the same as the length of cells.per.cluster")
    }
    
    cond.list <- list()
    cell.list <- list()
    for(i in seq_along(condition.props)){
        set.seed(42)
        block.cells <- cells.per.cluster[i]
        cell.list[[paste0("Block", i)]] <- block.cells
        
        block.cond <- rep("A", block.cells)
        block.a <- sample(1:block.cells, size=floor(block.cells*condition.props[i]))
        block.b <- setdiff(1:block.cells, block.a)
        block.cond[block.b] <- "B"
        cond.list[[paste0("Block", i)]] <- block.cond
    }
    blocks <- lapply(c(1:length(cell.list)), FUN=function(X) rep(paste0("B", X), cell.list[[X]]))
    rep.prop <- round(1/length(cells.per.cluster), 2)
    reps <- lapply(c(1:length(cell.list)), FUN=function(X) c(rep("R1", floor(cell.list[[X]] * rep.prop)),
                                                             rep("R2", floor(cell.list[[X]] * rep.prop)),
                                                             rep("R3", cell.list[[X]] - (2*floor(cell.list[[X]] * rep.prop)))))
    
    meta.df <- data.frame("Block"=unlist(blocks),
                          "Condition"=unlist(cond.list),
                          "Replicate"=unlist(reps))
    colnames(meta.df) <- c("Block", "Condition", "Replicate")
    rownames(meta.df) <- paste0("Cell", 1:nrow(meta.df))
    # define a "sample" as teh combination of condition and replicate
    meta.df$Sample <- paste(meta.df$Condition, meta.df$Replicate, sep="_")
    meta.df$Vertex <- c(1:nrow(meta.df))
    
    sim.sce <- SingleCellExperiment(assays=list(logcounts=sim.gex),
                                    reducedDims=list("PCA"=sim.pca$x))
    
    sim.mylo <- Milo(sim.sce)
    return(list("mylo"=sim.mylo, "meta"=meta.df))    
}
                     
message(paste0("Simulating ", opt$ncells, " across ", opt$clusters, " clusters."))
l.percluster <- rep(floor(opt$ncells/opt$clusters), opt$clusters)
n.fin <- opt$clusters - 1
l.percluster[opt$clusters] <- opt$ncells - sum(l.percluster[1:n.fin])

cond.props <- as.numeric(unlist(strsplit(opt$props, split=",", fixed=TRUE)))
sim.mylo <- simulate_discrete_data(n.clusters=opt$clusters, total.size=opt$ncells,
                                       cells.per.cluster=l.percluster, condition.props=cond.props)

ofile <- paste(opt$output, paste(paste0("Clusters", opt$clusters), paste0("Ncells", opt$ncells), opt$suffix, "discreteSim.RDS", sep="_"), sep="/")
message(paste0("Saving simulated data to: ", ofile))
saveRDS(sim.mylo, file=ofile)



