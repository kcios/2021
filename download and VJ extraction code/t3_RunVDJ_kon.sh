#!/bin/bash

#SBATCH --job-name=VDJ_Job
#SBATCH --time=72:00:00
#SBATCH --output=out.RunVDJ.%A_%a
#SBATCH --error=err.RunVDJ.%A_%a
#SBATCH --array=0-121
#SBATCH --mem=100250
#SBATCH --partition=rra
#SBATCH --qos=rra
#SBATCH --ntasks-per-node=20
#SBATCH --nodes=1

module load apps/python/3.8.5
IFS='#' read -a receptors_to_do < t3_rectodo.txt
#receptors_to_do=(TRA TRB TRD TRG IGH IGK IGL TRA_UM TRB_UM TRD_UM TRG_UM IGH_UM IGK_UM IGL_UM)

python3 -u "/work/pi_gblanck/Konrad/DLBCL_RNA/scripts/t3_findvdjum.py" ${receptors_to_do[$SLURM_ARRAY_TASK_ID]} '/work/pi_gblanck/Konrad/DLBCL_RNA/DLBCLRNAbams_Results/' '/work/pi_gblanck/Konrad/DLBCL_RNA/scripts/t3_vdjdb/'

echo 'done'