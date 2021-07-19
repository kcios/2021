#!/bin/bash
#SBATCH --job-name=TemplateSearchIgTcR
#SBATCH --time=24:00:00
#SBATCH --partition=rra
#SBATCH --qos=rra
#Notes below for how to estimate time based on what processes need to be run

# This entire job runs on 1 node, so this is effectively the cpu core limit.
#SBATCH --ntasks-per-node=24  
#SBATCH --mem=4096

#SBATCH --output=output.%j
#SBATCH --error=error.%j

#Max RAM usage ~7x the largest .tsv generated times number of cores; so 4 cores is 28*largest tsv

#160719 Use GNU parallel: Tange, O. (2020, September 22). GNU Parallel 20200922 ('Ginsburg'). Zenodo. https://doi.org/10.5281/zenodo.4045386

############################################################################################################################
#set file paths
bamFolder="/work/pi_gblanck/Konrad/CMI/CMI_bams"

#Set which IgTcrs to process
IgTcrs=("IGK" "IGH" "IGL" "TRA" "TRB" "TRD" "TRG" "IGH_UM" "IGK_UM" "IGL_UM" "TRA_UM" "TRB_UM" "TRD_UM" "TRG_UM")
#IgTcrs=("IGH" "IGK" "IGL" "TRA" "TRB" "TRD" "TRG" "IGH_UM" "IGK_UM" "IGL_UM" "TRA_UM" "TRB_UM" "TRD_UM" "TRG_UM")

#Set what needs to be executed, tsv2xlsx requires a lot of RAM so adjust accordingly
#imgtSearch is 20 minutes per file for all receptors combined
#extractUnmapped is about 8 minutes per file per core
#extractReceptors is about 2 minutes per file per core
#quickcheck is 1 second per file per core
#runSearches is 40 seconds per file per core for all receptors
#tsv2xlsx is 10 seconds per file per core
#so total time estimation is about 7 days for 5 TB/500 files
#Give extra 25% in case of hyperthreading vs physical cores
toDo=("quickcheck", "extractUnmapped", "extractReceptors","runSearches","tsv2xslx")
#todo=("imgtSearch")
#toDo=("quickcheck" "extractUnmapped" "extractReceptors" "runSearches" "tsv2xlsx" "imgtSearch" "emailStart" "emailEnd" "pigz" "copyToWos")

#don't set above the number of cores on the node, make sure if on circe, do NOT exceed ntasks-per-node above
#4 is best for scheduling
numCores="24"
#numCores=`nproc`

#comment/uncomment one of the lines below for cases where "chr" is required in the search script. It seems the downloaded bams have this chr prefix that the script will break without.
#chr=""
chr="chr"

email="kcios@usf.edu"
############################################################################################################################

module load apps/samtools/1.3.1
source /work/pi_gblanck/Konrad/CMI/scripts/t2_Module_Search_IgTcRFix.sh

