#!/bin/bash -l

#SBATCH --job-name=stage_1_offline_parallel
#SBATCH --nodes=1
#SBATCH --export=ALL
#SBATCH --output='slurm_output/slurm-%A_%a.out'


module load $6
cd $7
/usr/bin/time -v matlab -nodesktop -nodisplay -r "offline_parallel_stage_1($1,'$2','$5',$3,'$4','$8');exit;"
wait
echo "All Processess are Complete"
