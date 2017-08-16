#!/bin/bash
#SBATCH -J humanevol
#SBATCH -o /n/regal/debivort_lab/claire/humanevol/runevol.out
#SBATCH -e /n/regal/debivort_lab/claire/humanevol/runevol.err
#SBATCH -N 1 		
#SBATCH -c 32 		
#SBATCH -t 4-00:00 	
#SBATCH -p general 	
#SBATCH --mem=32000	
#SBATCH --mail-type=END
#SBATCH --mail-user=guerin.claire01@gmail.com

mkdir -p /scratch/$USER/$SLURM_JOB_ID

source new-modules.sh
module load python/2.7.6-fasrc01
srun -c 32 python runsim.py

rm -rf /scratch/$USER/$SLURM_JOB_ID