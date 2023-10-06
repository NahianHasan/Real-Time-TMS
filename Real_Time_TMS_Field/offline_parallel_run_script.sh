#!/bin/sh -l

#change the max wall time according to the specific cluster
stage_1_max_walltime=${12}
stage_2_max_walltime=${13}
stage_3_max_walltime=${14}
stage_4_max_walltime=${15}

#submit parallel jobs
while true; do
    DP=$(sbatch -A ${11} --cpus-per-task=$7 --time=$stage_1_max_walltime ${17}/offline_parallel_stage_1.sh $1 $2 $3 $5 $6 ${16} ${17} ${18})
    if [ "$?" = "0" ]; then
		break
	else
		sleep 600
	fi
done
while true; do
    CR=$(sbatch -A ${11} --dependency=afterany:${DP##* } --array=[1-$1] --cpus-per-task=$8 --time=$stage_2_max_walltime ${17}/offline_parallel_stage_2.sh $1 $3 $5 ${16} ${17})
    if [ "$?" = "0" ]; then
		break
	else
		sleep 600
	fi
done
while true; do
    CC=$(sbatch -A ${11} --dependency=afterany:${CR##* } --cpus-per-task=$9 --time=$stage_3_max_walltime ${17}/offline_parallel_stage_3.sh $1 $3 $5 ${16} ${17})
    if [ "$?" = "0" ]; then
		break
	else
		sleep 600
	fi
done
while true; do
    RC=$(sbatch -A ${11} --dependency=afterany:${CC##* } --cpus-per-task=${10} --time=$stage_4_max_walltime --array=[1-$1] ${17}/offline_parallel_stage_4.sh $1 $3 $5 ${16} ${17})
    if [ "$?" = "0" ]; then
		break
	else
		sleep 600
	fi
done
echo "Parallel jobs submitted"
