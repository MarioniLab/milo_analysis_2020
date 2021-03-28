### Run Milo vs MELD ###
suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(SingleCellExperiment)
  # library(scran)
  library(glue)
})

source('./benchmark_utils.R')

parser <- ArgumentParser()
parser$add_argument("data_RDS", type="character",
                    help = "path to RDS storing SingleCellExperiment object")
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
args <- parser$parse_args()

seed <- args$batch_seed
data_path <- args$data_RDS
k <- args$k
pop <- args$population
pop_enr <- args$pop_enrichment
be_sd <- args$batchEffect_sd

## Load data
print("Loading dataset...")
sce <- readRDS(data_path)

## Load coldata and PCA
outdir <- '/nfs/team205/ed6/data/milo_benchmark/synthetic_data/'
outprefix <- str_c("benchmark_embryo_pop_", pop, '_enr', pop_enr, "_seed", seed)
coldata <- read_csv(paste0(outdir, outprefix, ".coldata.csv")) %>% column_to_rownames()
X_pca <-read_csv(str_c(outdir, outprefix, "_batchEffect", be_sd, ".pca.csv")) %>% column_to_rownames()  

## Add reduced dim + coldata to sce
colData(sce) <- DataFrame(coldata)
reducedDim(sce, "pca_batch") <- as.matrix(X_pca)

## Run Milo and save
print("Running Milo...")
if (be_sd==0) {
  milo_res <- run_milo(sce, condition_col='synth_labels', sample_col="synth_samples",
                       reduced.dim = "pca_batch", d=30, k=k)
} else {
    ## Include batch in experimental design if batch effect is present
  milo_res <- run_milo(sce, condition_col='synth_labels', sample_col="synth_samples", batch_col="synth_batches",
                       reduced.dim = "pca_batch", d=30, k=k)
}

# extract results
milo_res$DAres[,"nhoodIndex"] <- unlist(nhoodIndex(milo_res$Milo))
milo_res$DAres[,'nh_size'] <- colSums(nhoods(milo_res$Milo))
milo_res$DAres <- annotateNhoods(milo_res$Milo, milo_res$DAres, coldata_col = "true_labels")
milo_outfile <- glue('/nfs/team205/ed6/data/milo_benchmark/MELDvsMilo/benchmark_embryo_pop_{pop}_enr{pop_enr}_seed{seed}_batchEffect{be_sd}.DAresults.milo.k{k}.csv')
print("Saving Milo output...")
write_csv(milo_res$DAres, milo_outfile)

## Run MELD and save
meld_call <- glue('python ./run_meld.py {seed} {pop} --pop_enrichment {pop_enr} --batchEffect_sd {be_sd} --k {k} --data_id embryo')
meld_outfile <- glue('/nfs/team205/ed6/data/milo_benchmark/benchmark_embryo_pop_{pop}_enr{pop_enr}_seed{seed}_batchEffect{be_sd}.DAresults.meld.csv')
meld_outfile_final <- glue('/nfs/team205/ed6/data/milo_benchmark/MELDvsMilo/benchmark_embryo_pop_{pop}_enr{pop_enr}_seed{seed}_batchEffect{be_sd}.DAresults.meld.k{k}.csv')
print("Running MELD...")
system(meld_call)
if (file.exists(meld_outfile)){
  rename_call <- glue("mv {meld_outfile} {meld_outfile_final}")
  print("Saving MELD file...")
  system(rename_call)
}


