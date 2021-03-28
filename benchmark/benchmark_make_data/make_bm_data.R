### Benchmarking methods with synthetic labels and batch effects on real data

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
parser$add_argument("batch_seed", type="integer",
                    help = "Seed 4 batch effect")
parser$add_argument("population", type="character",
                    help = "Which cell type is DA?")
parser$add_argument("--pop_enrichment", type="double", default=0.7,
                    help = "Max condition probability in DA population")
# parser$add_argument("--max_size", type="integer", default=1000,
#                     help = "Min number of cells in population to select")
parser$add_argument("--k", type="integer", default=50,
                    help = "K parameter")
parser$add_argument("--data_id", type="character", default="embryo",
                    help = "ID for the dataset used")
parser$add_argument("--make_batch_effect", type="character", default="yes",
                    help = "should synthetic batch effects be added? (yes/no)")
args <- parser$parse_args()

seed <- args$batch_seed
data_path <- args$data_RDS
k <- args$k
pop <- args$population
pop_enr <- args$pop_enrichment
data_id <- args$data_id
make_batch_effects <- args$make_batch_effect

## Load data
print("Loading dataset...")
sce <- readRDS(data_path)

## Select population to simulate DA by size and save
# sized_pops = names(table(sce$celltype))[table(sce$celltype) < max_size]
# pop = sample(sized_pops, 1)

if (str_detect(pop, "_")) {
  pop <- str_replace(pop, "_", " ")
}
sce <- add_synthetic_labels_pop(sce, pop=pop, pop_column = "celltype", seed=seed, pop_enr=pop_enr)
## Consider TRUE DA effect size +- 10%
if (pop_enr < 0.5) {
  da_lower <- pop_enr + (pop_enr/100)*10
  da_upper <- 1 - da_lower
} else {
  da_upper <- pop_enr - (pop_enr/100)*10
  da_lower <- 1 - da_upper
}
true_labels <- ifelse(sce$Condition2_prob < da_lower, "NegLFC", ifelse(sce$Condition2_prob > da_upper, "PosLFC", "NotDA"))
colData(sce)[["true_labels"]] <- true_labels
if (str_detect(pop, " ")) {
  pop <- str_replace(pop, " ", "_")
}

## Save coldata
outdir <-'/nfs/team205/ed6/data/milo_benchmark/synthetic_data/'
outprefix <- str_c("benchmark_", data_id, "_pop_", pop, '_enr', pop_enr, "_seed", seed)
coldata <- data.frame(colData(sce)) %>% rownames_to_column()
write_csv(coldata, str_c(outdir, outprefix, ".coldata.csv"))

## Simulate batch effects of different magnitude
set.seed(seed)
if (make_batch_effects=="yes") {
  print("Simulating batch effects...")
  bm_sce_ls <- lapply(c(0, 0.25, 0.5, 0.75, 1), function(sd){
    sce_be <- add_batch_effect(sce, batch_col = "synth_batches", norm_sd=sd)
    
    X_pca <- reducedDim(sce_be, "pca_batch")
    
    ## Save reduced dims
    write_csv(as.data.frame(X_pca) %>% rownames_to_column(), str_c(outdir, outprefix, "_batchEffect", sd, ".pca.csv"))
  })
} else {
  X_pca <- reducedDim(sce, "pca.corrected")
  ## Save reduced dims
  write_csv(as.data.frame(X_pca) %>% rownames_to_column(), str_c(outdir, outprefix, "_batchEffect0.pca.csv"))
}
