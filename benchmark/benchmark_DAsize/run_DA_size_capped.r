### Run DA methods in R ###

suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(SingleCellExperiment)
  library(scran)
})

source('./benchmark_utils.R')

parser <- ArgumentParser()
parser$add_argument("data_RDS", type="character",
                    help = "path to RDS storing SingleCellExperiment object")
parser$add_argument("method", type="character",
                    help = "DA method to use")
parser$add_argument("batch_seed", type="integer",
                    help = "Seed 4 batch effect")
parser$add_argument("population", type="character",
                    help = "Which cell type is DA?")
parser$add_argument("--pop_enrichment", type="double", default=0.7,
                    help = "Max condition probability in DA population")
parser$add_argument("--batchEffect_sd", type="double", default=0,
                    help = "SD of batch effect")
parser$add_argument("--k", type="integer", default=50,
                    help = "K parameter")
parser$add_argument("--data_id", type="character", default="embryo",
                    help = "ID for the dataset used")
parser$add_argument("--MNN_correct", type="character", default="no",
                    help = "should I use MNN corrected version")
args <- parser$parse_args()

seed <- args$batch_seed
data_path <- args$data_RDS
k <- args$k
pop <- args$population
pop_enr <- args$pop_enrichment
be_sd <- args$batchEffect_sd
DA_method <- args$method
data_id <- args$data_id
mnn_correct <- args$MNN_correct

## Load data
print("Loading dataset...")
sce <- readRDS(data_path)

## Load coldata and PCA
outdir <- '/nfs/team205/ed6/data/milo_benchmark/synthetic_data/'
outprefix <- str_c("benchmark_", data_id, "_pop_", pop, '_enr', pop_enr, "_seed", seed, "_capped")
coldata <- read_csv(paste0(outdir, outprefix, ".coldata.csv")) %>% column_to_rownames()
if (mnn_correct=="no") {
  X_pca <-read_csv(str_c(outdir, outprefix, "_batchEffect", be_sd, ".pca.csv")) %>% column_to_rownames()  
} else {
  X_pca <-read_csv(str_c(outdir, outprefix, "_batchEffect", be_sd, ".MNNcorrected.pca.csv")) %>% column_to_rownames()
}


## Find DA probability x cell
tol_dataset <- list(cluster=0.8, branching=6, linear=6, 
                    cluster3x=0.5, branching3x=6, linear3x=6,
                    embryo=1) ## cydar radius picked w heuristic
bm_params = list(
  milo = list(k=k),
  milo_batch = list(k=k),
  meld = list(k=k),
  daseq = list(k.vec=seq(k, 500, 50)),
  louvain = list(k=k),
  louvain_batch = list(k=k),
  cydar = list(tol=tol_dataset[[data_id]], downsample=3),
  cydar_batch = list(tol=tol_dataset[[data_id]], downsample=3)
  )

## Run DA method
out <- runDA(sce, X_pca, coldata = coldata, method = DA_method, params = bm_params)

## Save output 
bm_outdir <-'/nfs/team205/ed6/data/milo_benchmark/'
if (mnn_correct=="no") {
  write_csv(out, str_c(bm_outdir, outprefix, "_batchEffect", be_sd, ".DAresults.", DA_method, ".csv"))
} else {
  write_csv(out, str_c(bm_outdir, outprefix, "_batchEffect", be_sd, ".MNNcorrected.DAresults.", DA_method, ".csv"))
  } 
