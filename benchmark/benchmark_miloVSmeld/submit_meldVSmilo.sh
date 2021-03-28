#!/bin/bash

## Activate conda env
conda activate milo_bm

outdir=/nfs/team205/ed6/bin/milo_benchmarking/outfiles
data_id=embryo
pop_file=/nfs/team205/ed6/data/milo_benchmark/pop_sample_2_clean.txt
pops=$(cat $pop_file | while read pop; do echo $pop; done)
seed=43

## Run
for pop in $pops
        do
          for enr in 0.7 0.8 0.9
            do
            for be_sd in 0 0.25 0.5
                do
                        echo "Rscript ./run_miloVSmeld.R --pop_enrichment $enr --batchEffect_sd $be_sd --k 50 /nfs/team205/ed6/data/milo_benchmark/embryo_data_bm.RDS $seed $pop" | bsub -o ${outdir}/miloVSmeld_bm_data_${pop}_${seed}.${enr}.out -e ${outdir}/miloVSmeld_bm_data_${pop}_${seed}.${enr}.err -G team283 -R"select[mem>3500] rusage[mem=3500]" -M3500
        done
    done
    done
