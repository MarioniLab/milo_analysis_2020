### Scalability analysis downsampling liver data ##
suppressPackageStartupMessages({
  library(tidyverse)
  library(irlba)
  library(DropletUtils)
  library(scater)
  library(scran)
  library(SingleCellExperiment)
  library(miloR)
  library(parallel)
  library(pbmcapply)
})

liver_milo <- readRDS("/nfs/team205/ed6/data/Ramachandran2019_liver/tissue_milo.RDS")
liver_milo <- Milo(as(liver_milo, "SingleCellExperiment"))

## Downsample milo obj
downsample_milo <- function(liver_milo, prop=0.1, seed=42){
  set.seed(seed)
  smp_cells <- sample(colnames(liver_milo), size = round(ncol(liver_milo)*prop))
  dsmp_milo <- liver_milo[,smp_cells]
  return(dsmp_milo)
}

seeds <- 42 + seq(1,3)
props <- seq(0.05,0.95, by=0.05)

param_mat <- cbind(unlist(lapply(props, rep, 3)), rep(seeds, length(props)))

dsmp_milo_ls <- lapply(1:nrow(param_mat), function(x) downsample_milo(liver_milo, prop=param_mat[x,1], seed=param_mat[x,2]))

#### Make design matrix
liver_meta <- as.tibble(colData(liver_milo)[,c("dataset","condition")]) 
liver_meta <- distinct(liver_meta) %>%
  mutate(condition=factor(condition, levels=c("Uninjured", "Cirrhotic"))) %>%
  column_to_rownames("dataset")

## Run Milo analysis
run_milo <- function(milo, k = 30){
  milo <- buildGraph(milo, d = 11, k=k)
  milo <- makeNhoods(milo, k=k, d=11, refined=TRUE)
  milo <- countCells(milo, 
                     samples = "dataset", 
                     meta.data = data.frame(colData(milo)[,c("dataset","condition")]) 
                     )
  milo_res <- testNhoods(milo, 
                         design = ~ condition, 
                         design.df = liver_meta, 
                         fdr.weighting = "k-distance")
  return(milo_res)
}

measure_milo <- function(milo){
  run_time <- system.time(res <- run_milo(milo))["elapsed"]  
  n_nhoods <- nrow(res)
  return(list(n_nhoods=n_nhoods, run_time=run_time))
  }


results <- pbmclapply(1:length(dsmp_milo_ls), 
                    function(i){
                      print(paste0("Downsample no.", i))
                      measure_milo(dsmp_milo_ls[[i]])
                    }, mc.cores=10, mc.preschedule = FALSE)


saveRDS(results, "/nfs/team205/ed6/data/Ramachandran2019_liver/scalability_results.RDS")
results <- readRDS( "/nfs/team205/ed6/data/Ramachandran2019_liver/scalability_results.RDS")

results
results_df <- lapply(seq_along(results), function(i) data.frame(results[[i]])) %>%
  purrr::reduce(bind_rows) %>%
  mutate(prop = param_mat[1:50,1], seed = param_mat[1:50,2]) 

results_df %>% 
  mutate(n_cells = round(ncol(liver_milo)*prop)) %>%
  ggplot(aes(n_cells, run_time, fill=n_nhoods)) + 
  geom_point(size=3, color="black", pch=21) +
  scale_fill_viridis_c(name="# nhoods") +
  xlab("# cells") + ylab("Run time (s)") +
  theme_clean(base_size = 18) +
  ggtitle("Downsampled liver dataset") +
  ggsave("/nfs/team205/ed6/data/Ramachandran2019_liver/scalability_results.png", width=8, height = 7)

write_csv(results_df, "/nfs/team205/ed6/data/Ramachandran2019_liver/scalability_results.csv")


