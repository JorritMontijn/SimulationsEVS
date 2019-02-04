#!/bin/sh

#SBATCH --job-name="IndInpV1V2Ori5Noise5"
#SBATCH --partition=mono-shared
#SBATCH --time=0-12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000

#SBATCH --output=/home/montijn/IndInpV1V2Ori5Noise5.out
#SBATCH --error=/home/montijn/IndInpV1V2Ori5Noise5.err


module load foss/2016a matlab/2016b
srun ./DynSimCompiledIndInp time=0-12:00:00,conn=sConnSimil2_Ret32Col60N2400S725760_2018-11-23.mat,stim=sStim2_SquareGratingExpRet32Noise5Ori5Drift2_x2R1_2018-08-08.mat,idx=$1,att=0_2,tag=IndInpV1V2Ori5Noise5
