## Benchmarking scripts

All analyses and steps are described in `Fig2_benchmark_results.Rmd`

### Contents

* `benchmark_make_data/`: scripts to construct synthetic labels on existing datasets
* `benchmark_main/`: scripts to run main benchmarking analysis (4 methods, 4 datasets, + batch effects)
* `benchmark_DAsize`: scripts to run benchmark on DA regions of controlled size (Supp Fig 6)
* `benchmark_Kselection`: scripts to run milo on semi-synthetic dataset controlling K parameter (Supp Notes Fig 1)
* `benchmark_MNN`: scripts to run benchmark on data with MNN corrected synthetic batch effects (Supp Fig 10)
* `benchmark_miloVSmeld`: scripts to run comparison between milo and MELD (Supp Fig 7)
* `benchmark_utils.R`: utility functions to generate synthetic condition labels and to run DA methods for benchmarking
* `run_DA_R.r`: utility script to run DA methods on pre-computed synthetic labels and datasets