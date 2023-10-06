#!/bin/bash
# FILENAME:  TMS_simulation_ground_truth

#SBATCH --nodes=1
#SBATCH --export=ALL
#SBATCH --job-name 'TMS_grnd_trth'
#SBATCH --output='./slurm_output/slurm-%A_%a.out'
#SBATCH --error='./slurm_output/slurm-%A_%a.out'

# Loads Matlab and sets the application up
module load ${10}

ulimit -c 0                    # no core dumps
ulimit -s unlimited            # unlimited stack

cd $1
matlab -nodisplay -nodesktop -r "ground_truth('$1',${SLURM_ARRAY_TASK_ID},'$2','$3','$4','$5','$6',$7,'$8','$9'); exit;"&
wait

echo "All Processess are Complete"


