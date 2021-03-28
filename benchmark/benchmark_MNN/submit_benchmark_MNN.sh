#!/bin/bash

outdir=/nfs/team205/ed6/bin/milo_benchmarking/outfiles
## Activate conda env
conda activate milo_bm

data_id=embryo


pop_file=/nfs/team205/ed6/data/milo_benchmark/pop_sample_2_clean.txt
pops=$(cat $pop_file | while read pop; do echo $pop; done)
R_methods=$(for m in milo daseq cydar louvain milo_batch louvain_batch cydar_batch; do echo $m; done)
batch_vec=$(for m in 0.5 0.75 1; do echo $m; done)
k=50
pop_enr=0.8
seed=43
output_dir=/nfs/team205/ed6/data/milo_benchmark/

## Run
for pop in $pops
        do
        for batch_sd in $batch_vec
            do
            for method in $R_methods
            do
            echo "Rscript ./run_DA_R.r /nfs/team205/ed6/data/milo_benchmark/${data_id}_data_bm.RDS $method $seed $pop --pop_enrichment $pop_enr --data_id $data_id --k $k --batchEffect_sd $batch_sd --MNN_correct yes" | \
                                bsub -G team283  -o ${outdir}/milo_bm_MNN_${data_id}_${seed}_${method}_${pop}_${batch_sd}.out -e ${outdir}/milo_bm_MNN_${data_id}_${seed}_${method}_${pop}_${batch_sd}.err -R"select[mem>3500] rusage[mem=3500]" -M3500 
                    done
                done
            done
        done
done
    
