#!/bin/sh

#SBATCH --job-name="Ori5Noise02"
#SBATCH --partition=mono-shared
#SBATCH --time=0-12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000

#SBATCH --output=/home/montijn/Ori5Noise02.out
#SBATCH --error=/home/montijn/Ori5Noise02.err


module load foss/2016a matlab/2018b
srun ./SimCompiled time=0-12:00:00,conn=Conn256N1200_2020-07-16.mat,stim=Ret256Noise0.2Ori5_x2R1_2020-07-17.mat,idx=$1,tag=Ori5Noise02
