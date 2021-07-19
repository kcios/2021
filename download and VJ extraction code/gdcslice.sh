#!/bin/bash
#SBATCH --job-name=GDC-Download-sliced
#SBATCH --time=1:00:00
#SBATCH --ntasks=4
#SBATCH --partition=rra
#SBATCH --qos=rra
#SBATCH --mem-per-cpu=2048
#SBATCH --output=output.%j

# channge this number to the total number of files in your manifest, ie the number 372 do not change the %100
#SBATCH --array=0-200%100  # arent bash arrays 0 indexed? possible bug? dont think it matters though.


module purge
module add apps/gdc-client/1.3.0 # is this even needed? Dont think so, but will keep it for fun.

#path to manifest, change this to your manifest location
PathToManifest="/work/pi_gblanck/Konrad/DLBCL_RNA/manifest.txt"

#path to download token, change this to your token location
Token="/work/pi_gblanck/Konrad/DLBCL_RNA/gdctoken.txt"

# change this path tp the path to you want the files to go to
OutputFolder="/work/pi_gblanck/Konrad/DLBCL_RNA/DLBCLRNAbams"

mapfile -t myArray < $PathToManifest

NumberOfBams=${#myArray[@]}
echo $NumberOfBams
InputString=${myArray[$SLURM_ARRAY_TASK_ID]}
ID=$(cut -d' ' -f1 <<< $InputString)
NAME=$(cut -d' ' -f2 <<< $InputString)

APItext="https://api.gdc.cancer.gov/slicing/view/$ID?region=chr14:21521904-22652132&region=chr7:142289011-142913287&region=chr7:38240024-39368055&region=chr14:105486437-106879844&region=chr2:88857361-90238368&region=chr22:22026076-22922913&region=chr6:29844528-33100696&region=unmapped"
#                                                             TRA / TRD                       TRB                              TRG                          IGH                                 IGK                               IGL                            HLA
token=$(<$Token)
curl --header "X-Auth-Token: $token" $APItext --output $OutputFolder/sliced_$NAME

#token=$(<$Token)
#curl --header "X-Auth-Token: $token" 'https://api.gdc.cancer.gov/slicing/view/df80679e-c4d3-487b-934c-fcc782e5d46e?region=chr14:21621904-22552132&region=chr3:10000-20000' --output get_regions_slice.bam

