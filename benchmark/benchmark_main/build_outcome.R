## Compute benchmark outcomes ##  
suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(SingleCellExperiment)
  library(scran)
})

source('./benchmark_utils.R')

## Calculate outcome by setting a probability threshold for true DA (da_upper)
.outcome_by_prob <- function(benchmark_df, da_upper){
  DA_thresh <- 1 - da_upper
  benchmark_df <- benchmark_df %>%
    mutate(true = ifelse(true_prob < DA_thresh, "NegLFC", "NotDA")) 
  if (benchmark_df$method[1]=="meld") {
    benchmark_df <- mutate(benchmark_df, pred = ifelse(Condition2 < DA_thresh, "NegLFC", "NotDA")) 
  }
  ## Check if conditions were swapped in test
  pred_cor <- benchmark_df %>%
    mutate(pred=factor(pred, levels=c("NegLFC", "NotDA", "PosLFC"), ordered = TRUE)) %>%
    mutate(pred=as.numeric(pred)) %>%
    summarise(cor(pred, true_prob))
  ## Swap outcomes if so
  if (!is.na(pred_cor)) {
    if (pred_cor < -0.1) {
      benchmark_df <- mutate(benchmark_df, pred = ifelse(pred=="NegLFC", "PosLFC", ifelse(pred=="PosLFC", "NegLFC", "NotDA")))
    }
  }
  calculate_outcome(benchmark_df)
}

## Read outputs
outdir <- "/nfs/team205/ed6/data/milo_benchmark/embryo/"
res_files <- list.files(outdir, pattern=".DAresults")
res_files_full <- list.files(outdir, pattern=".DAresults", full.names = TRUE)
# in_files <- list.files("~/data/milo_benchmark/synthetic_data/", pattern="coldata.csv")

## Make data frame w benchmark parameters
out_meta_df <- data.frame(file_id = str_remove_all(res_files, "benchmark_embryo_pop_|.csv")) %>%
  separate(col = file_id, sep = ".DAresults.", into=c("file_id", "method")) %>%
  separate(col = file_id, sep = "_enr", into=c("pop", "file_id")) %>%
  separate(col = file_id, sep = "_", into=c("enr", "seed", "batchEffect")) %>%
  mutate(seed=str_remove(seed, "seed"), batchEffect=str_remove_all(batchEffect, "batchEffect")) 

## Read benchmark results w/o batch effects
# using different probability thresholds
no_batch_files <- str_detect(res_files_full, "batchEffect0.DAresult") & !str_detect(res_files_full, "louvain_batch|cydar_batch|milo_batch")
prob_thresh_vec <- seq(0.5, 0.9, 0.05)
outcome_df <- lapply(seq_along(res_files_full[no_batch_files]), function(i){ 
  print(paste("Outcome no. ", i))
  benchmark_df <- read_csv(res_files_full[no_batch_files][i])
  pop_enr <- as.numeric(out_meta_df[no_batch_files,][i,"enr"])
  lapply(prob_thresh_vec[prob_thresh_vec < pop_enr], function(x){
    benchmark_df %>%
      .outcome_by_prob(da_upper = x) %>%
      mutate(DA_thresh=x) %>%
      ungroup() %>%
      dplyr::select(- method) %>%
      bind_cols(out_meta_df[no_batch_files,][i,]) %>%
      filter(DA_thresh < as.numeric(enr))
  }) %>%
    purrr::reduce(bind_rows)
}) %>%
  purrr::reduce(bind_rows)

write_csv(outcome_df, "/nfs/team205/ed6/data/milo_benchmark/outcome_embryo_v3_noBatch.csv")

# Read benchmark results w batch effects
# using probability threshold = 0.6
outcome_df <- lapply(seq_along(res_files_full[!no_batch_files]), function(i){ 
  print(paste("Outcome no. ", i))
  benchmark_df <- read_csv(res_files_full[!no_batch_files][i])
  pop_enr <- as.numeric(out_meta_df[!no_batch_files,][i,"enr"])
  benchmark_df %>%
      .outcome_by_prob(da_upper = 0.6) %>%
      mutate(DA_thresh=0.6) %>%
      ungroup() %>%
      dplyr::select(- method) %>%
      bind_cols(out_meta_df[!no_batch_files,][i,]) %>%
      filter(DA_thresh < as.numeric(enr))
  }) %>%
  purrr::reduce(bind_rows)

write_csv(outcome_df, "/nfs/team205/ed6/data/milo_benchmark/outcome_embryo_v3_wBatch.csv")

