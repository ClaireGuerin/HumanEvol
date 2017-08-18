#!/bin/bash
#SBATCH -J humanevol
#SBATCH -o /n/regal/debivort_lab/claire/humanevol/runevol.out
#SBATCH -e /n/regal/debivort_lab/claire/humanevol/runevol.err
#SBATCH -N 1 		
#SBATCH -c 32 		
#SBATCH -t 4-00:00 	
#SBATCH -p serial_requeue	
#SBATCH --mem=32000
#SBATCH --mail-type=END
#SBATCH --mail-user=guerin.claire01@gmail.com

source new-modules.sh
module load python/2.7.6-fasrc01
	
while IFS='' read -r line || [[ -n $line ]]; do
    #echo $line
	mkdir -p /scratch/$USER/$SLURM_JOB_ID
	srun -c 32 python -u runsim.py $line
	sleep 1 # pause to be kind to the scheduler
	rm -rf /scratch/$USER/$SLURM_JOB_ID
done < "$1"
