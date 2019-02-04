#!/bin/sh

#SBATCH --job-name="C2ExpRet32Col180Ori2Noise063"
#SBATCH --partition=mono-shared
#SBATCH --time=0-12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000

#SBATCH --output=/home/montijn/C2ExpRet32Col180Ori2Noise063.out
#SBATCH --error=/home/montijn/C2ExpRet32Col180Ori2Noise063.err


module load foss/2016a matlab/2016b
srun ./DynSimCompiled time=0-12:00:00,conn=sConnSimil2_Ret32Col180N4800S1397760_2018-04-09.mat,stim=sStim2_SquareGratingExpRet32Noise0.63Ori2Drift2_x2R1_2018-06-27.mat,idx=$1,att=0_2,tag=C2ExpRet32Col180Ori2Noise063
