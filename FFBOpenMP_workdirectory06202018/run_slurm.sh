#!/bin/bash
#SBATCH --job-name="ffb0"
#SBATCH --output="ffb0.%j.%N.out"
#SBATCH --partition=compute
#SBATCH --export=ALL
#SBATCH --nodes 1
#SBATCH --ntasks-per-node=24
#SBATCH -t 10:30:00
#SBATCH -A cno104 
#SBATCH --mail-user=alexander.alekseenko@csun.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end


export OMP_NUM_THREADS=24 

cd $SLURM_SUBMIT_DIR

./ffbM650.a >outM650ffbM41tr3_TREF10000_Decomp_00.txt
