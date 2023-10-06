#!/bin/bash -l

#SBATCH --job-name=stage_3_offline_parallel
#SBATCH --nodes=1
#SBATCH --export=ALL
#SBATCH --output='slurm_output/slurm-%A_%a.out'


module load $4
cd $5
/usr/bin/time -v matlab -nodesktop -nodisplay -r "offline_parallel_stage_3($1,$2,'$3');exit;"
wait
echo "All Processess are Complete"
