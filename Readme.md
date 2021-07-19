## Code for DLBCL project


1. Run sbatch gdcslice.sh
	a. This will download slices from your manifest. 
	b. Get manifest and token file beforehand as .txt.
	c. Create a bams folder to hold the files, set the script to point to it.
2. Run sbatch t2_Module_Search_IgTcR_header.sh
	a. First step in processing the bams, will create a ton of files.
	b. Edit the file paths of this script at the top (bam folder) and bottom (script folder).
	c. Make sure GNU parallel is installed.
	d. Copy all t2 files over.
	e. Takes like 1 hour for 800ish bam files, of DNA and RNA types.
3. Run t3_set_task_items.py
	a. Edit the path at the top of the python file to match the results folder generated from the step before.
	b. Note the array setting printed out to put into the slurm config file, t3_Run_VDJ_kon.sh
4. Run sbatch t3_Run_VDJ_kon.sh
	a. Will pull stuff out of bams into csvs, also will find the matching V/J/D etc sequence.
	b. Need to edit the file paths in this script. Rest of them should be good.
	c. Remember to place the vdjdb folder into the proper directory you indicate.
5. t4_pre.py to convert csv to xlsx
	a. Don’t have to actually use the xlsx files, csv are provided as well. Sometimes the data will be too large for xlsx anyways.
	b. Will combine the partitioned files from the previous analysis, based on receptor type.
6. sbatch t4_start.sh shell script
	a. Will compile all the data and filter it.
	b. Edit cancer in VDJrecord, filepaths in t4 start
	c. Download sample.tsv, and put it in the results from t4_pre folder. It is under the biospecimin download button on GDC.
Can run Physiochem or not.