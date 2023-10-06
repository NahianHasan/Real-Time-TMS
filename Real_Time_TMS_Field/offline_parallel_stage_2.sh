#!/bin/bash -l


#SBATCH --job-name=stage_2_offline_parallel
#SBATCH --nodes=1
#SBATCH --export=ALL
#SBATCH --output='slurm_output/slurm-%A_%a.out'


module load $4
cd $5
/usr/bin/time -v matlab -nodesktop -nodisplay -r "offline_parallel_stage_2($1,$2,${SLURM_ARRAY_TASK_ID},'$3');exit;"
wait
echo "All Processess are Complete"
