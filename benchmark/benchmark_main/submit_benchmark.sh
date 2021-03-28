#!/bin/bash

outdir=/nfs/team205/ed6/bin/milo_benchmarking/outfiles
## Activate conda env
conda activate milo_bm

data_id=$1

if [[ "$data_id" == "cluster" ]]
    then 
    pops=$(seq 1 1 3)
    R_methods=$(for m in milo daseq cydar louvain; do echo $m; done)
    batch_vec=0
    k=30
elif [[ "$data_id" == "linear" ]]
    then 
    pops=$(for p in $(seq 1 1 7); do echo M$p; done)
    R_methods=$(for m in milo daseq cydar louvain; do echo $m; done)
    batch_vec=0
    k=30
elif [[ "$data_id" == "branching" ]]
    then 
    pops=$(for p in $(seq 1 1 10); do echo M$p; done)
    R_methods=$(for m in milo daseq cydar louvain; do echo $m; done)
    batch_vec=0
    k=30
elif [[ "$data_id" == "cluster3x" ]]
    then 
    pops=$(for p in $(seq 1 1 3); do echo B$p; done)
    R_methods=$(for m in milo daseq cydar louvain; do echo $m; done)
    batch_vec=0
    k=30
elif [[ "$data_id" == "linear3x" ]]
    then 
    pops=$(for p in $(seq 1 1 7); do echo M$p; done)
    R_methods=$(for m in milo daseq cydar louvain; do echo $m; done)
    batch_vec=0
    k=30
elif [[ "$data_id" == "branching3x" ]]
    then 
    pops=$(for p in $(seq 1 1 10); do echo M$p; done)
    R_methods=$(for m in milo daseq cydar louvain; do echo $m; done)
    batch_vec=0
    k=30
elif [[ "$data_id" == "embryo" ]]
    then 
    pop_file=/nfs/team205/ed6/data/milo_benchmark/pop_sample_2_clean.txt
    pops=$(cat $pop_file | while read pop; do echo $pop; done)
    R_methods=$(for m in milo daseq cydar louvain milo_batch louvain_batch cydar_batch; do echo $m; done)
    batch_vec=$(for m in 0 0.25 0.5 0.75 1; do echo $m; done)
    k=50
fi

output_dir=/nfs/team205/ed6/data/milo_benchmark/

## Run
for pop in $pops
        do
        for pop_enr in $(seq 0.75 0.1 0.95)
            do
            for seed in $(seq 43 1 45)
                do
                for batch_sd in $batch_vec
                    do
#                         ## Check if MELD file exists
#                         if [ -f $output_dir/benchmark_${data_id}_pop_${pop}_enr${pop_enr}_seed${seed}_batchEffect${batch_sd}.DAresults.meld.csv ]; then
#                                 echo "Output file exists"
#                             else
#                                 echo "python ./run_meld.py $seed $pop --pop_enrichment $pop_enr --batchEffect_sd $batch_sd --k $k --data_id $data_id" | bsub -G team283 -o ${outdir}/milo_bm_${data_id}_${seed}_meld_${pop}_${batch_sd}.out -e ${outdir}/milo_bm_${data_id}_${seed}_meld_${pop}_${batch_sd}.err -R"select[mem>3500] rusage[mem=3500]" -M3500
#                             fi
                        for method in $R_methods
                            do
                            # Check if outfile exists already
                            if [ -f $output_dir/benchmark_${data_id}_pop_${pop}_enr${pop_enr}_seed${seed}_batchEffect${batch_sd}.DAresults.${method}.csv ]; then
                                echo "Output file exists"
                            else
                                echo "Rscript ./run_DA_R.r /nfs/team205/ed6/data/milo_benchmark/${data_id}_data_bm.RDS $method $seed $pop --pop_enrichment $pop_enr --data_id $data_id --k $k --batchEffect_sd $batch_sd" | \
                                bsub -G team283  -o ${outdir}/milo_bm_${data_id}_${seed}_${method}_${pop}_${batch_sd}.out -e ${outdir}/milo_bm_${data_id}_${seed}_${method}_${pop}_${batch_sd}.err -R"select[mem>3500] rusage[mem=3500]" -M3500 
                            fi
                    done
                done
            done
        done
done
    
