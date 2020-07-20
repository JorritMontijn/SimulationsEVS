#!/bin/sh

#SBATCH --job-name="Ori5Noise20"
#SBATCH --partition=mono-shared-EL7
#SBATCH --time=0-12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10000
#SBATCH --mail-user=Jorrit.Montijn@unige.ch
#SBATCH --mail-type=FAIL
#SBATCH --output=/home/users/m/montijn/SimulationsEVS/Ori5Noise20.out
#SBATCH --error=/home/users/m/montijn/SimulationsEVS/Ori5Noise20.err



module load foss/2016a matlab/2018b
srun ./SimCompiled time=0-12:00:00,conn=Conn256N1200_2020-07-16.mat,stim=Ret256Noise2.0Ori5_x2R1_2020-07-17.mat,idx=$1,tag=Ori5Noise20
