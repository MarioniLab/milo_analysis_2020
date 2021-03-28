#!/bin/bash

outdir=/nfs/team205/ed6/bin/milo_benchmarking/outfiles
## Activate conda env
conda activate milo_bm

data_id=embryo
# pop_file=/nfs/team205/ed6/data/milo_benchmark/pop_sample_2_clean.txt
# pops=$(cat $pop_file | while read pop; do echo $pop; done)
R_methods=$(for m in milo daseq cydar louvain; do echo $m; done)
pop_enr=$(for m in 0.75 0.85 0.95; do echo $m; done)
k=50
seed=$(for m in 43 44 45; do echo $m; done)
output_dir=/nfs/team205/ed6/data/milo_benchmark/

## Run
for pop in Erythroid2 Caudal_neurectoderm
        do
        for s in $seed
            do
            for method in $R_methods
                do
                for enr in $pop_enr
                do
                for a_logit in 03 05 07
                do
            echo "Rscript ./run_DA_R.r /nfs/team205/ed6/data/milo_benchmark/${data_id}_data_bm.RDS $method $s $pop --pop_enrichment $enr --data_id $data_id --k $k --capped yes --logit_param $a_logit" | \
                                bsub -G team283  -o ${outdir}/milo_bm_capped_${data_id}_${s}_${method}_${pop}.out -e ${outdir}/milo_bm_capped_${data_id}_${s}_${method}_${pop}.err -R"select[mem>3500] rusage[mem=3500]" -M3500 
                    done
                done
            done
        done
done
    
