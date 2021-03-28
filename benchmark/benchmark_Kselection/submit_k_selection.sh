#!/bin/bash

## Activate conda env
conda activate milo_bm

outdir=/nfs/team205/ed6/bin/milo_benchmarking/outfiles
data_id=embryo
pop_file=/nfs/team205/ed6/data/milo_benchmark/pop_sample_2_clean.txt
pops=$(cat $pop_file | while read pop; do echo $pop; done)
k_vec=$(seq 10 10 100)

output_dir=/nfs/team205/ed6/data/milo_benchmark/

## Run
for pop in $pops
        do
          for k in $k_vec
            do
            if [ -f $output_dir/benchmark_${data_id}_pop_${pop}_enr0.8_seed43_batchEffect0.DAresults.milo.k${k}.csv ]; then
                                echo "Output file exists"
                            else
                        echo "Rscript ../run_DA_R.r /nfs/team205/ed6/data/milo_benchmark/${data_id}_data_bm.RDS milo 43 $pop --pop_enrichment 0.8 --data_id $data_id --k $k --batchEffect_sd 0; mv $output_dir/benchmark_${data_id}_pop_${pop}_enr0.8_seed43_batchEffect0.DAresults.milo.csv $output_dir/benchmark_${data_id}_pop_${pop}_enr0.8_seed43_batchEffect0.DAresults.milo.k${k}.csv " | \
                        bsub -G team283  -o ${outdir}/milo_bm_${data_id}_43_milo_k${k}.out -e ${outdir}/milo_bm_${data_id}_43_milo_k${k}.err -R"select[mem>3500] rusage[mem=3500]" -M3500 
                        fi
        done
    done

