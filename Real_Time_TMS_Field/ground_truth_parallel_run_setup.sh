#!/bin/sh -l

#change the max wall time according to the specific cluster
max_walltime=${13}


while true; do
    DP=$(sbatch -A ${11} --array=[1-${10}] --cpus-per-task=${12} --time=$max_walltime $1/ground_truth_cluster_run.sh $1 $2 $3 $4 $5 $6 $7 $8 $9 ${14})
    if [ "$?" = "0" ]; then
	    break
    else
	    sleep 600
    fi
done

echo "Parallel jobs submitted"


