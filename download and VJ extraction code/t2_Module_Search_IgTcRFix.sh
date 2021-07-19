imgtSearch="/work/pi_gblanck/old_blanck_group/shared/imgtSearchIgTcr.php"
searchScript="/work/pi_gblanck/Konrad/CMI/scripts/t2_SearchReadsFix.sh"

#Set bin paths
tsv2xlsx="/work/pi_gblanck/old_blanck_group/shared/tsv2xlsx.py"
#parallelBin="/home/t/thuda/bin" #notused bc installed personally
#parallel --citation
phpBin="/work/pi_gblanck/old_blanck_group/shared/bin/php"
pigz="/work/pi_gblanck/old_blanck_group/shared/bin/pigz"

# this is a sort of manager script that takes a bunch of tasks and then executes them.

# no configuration below this line !!!

module purge
module add apps/python/2.7.11
module load apps/samtools/1.3.1

samplesName="${bamFolder##*/}"

#Start the processing
echo "Starting up..."
echo `date +"%a %x %X"`

#Make results folder
resultsDir="${bamFolder}_Results"
echo bamFolder is ${bamFolder}
#resultsDir="${HOME}/${samplesName}_Results"

##############################################################

if [[ "${toDo[*]}" == *emailStart* ]]; then
	date +"%a %x %X" | mail -s "Script started" $email
fi

##############################################################

#Make the results folders for each receptor
echo "Making folders..."
echo `date +"%a %x %X"`

mkdir -p $resultsDir

for IgTcr in "${IgTcrs[@]}"
do
	resultsDirIgTcr="${resultsDir}/$IgTcr"

	mkdir -p $resultsDirIgTcr
done

##############################################################

echo "Fetching file list..."
echo `date +"%a %x %X"`

#Get list of file paths
#Check if Ig/TCR array requires unmapped reads
if [[ "${toDo[*]}" == *extractUnmapped* ]]; then
	if [[ "${IgTcrs[*]}" == *UM* ]]; then
		numBams=`find $bamFolder -type f -name "*.bam" ! -name "*unmapped*"  ! -name "*TRA.bam" ! -name "*TRB.bam" ! -name "*TRD.bam" ! -name "*TRG.bam" ! -name "*IGH.bam" ! -name "*IGK.bam" ! -name "*IGL.bam" -printf '.' | wc -c`
		numUM=`find $bamFolder -type f -name "*unmapped.bam" -printf '.' | wc -c`
		echo "Found $numBams bam files and $numUM unmapped bam files..."

		if [ "$numBams" -ne "$numUM" ]; then
			echo "Extracting unmapped regions..."
			fileListM=`find ${bamFolder} -type f -name "*.bam" ! -name "*unmapped*"  ! -name "*TRA.bam" ! -name "*TRB.bam" ! -name "*TRD.bam" ! -name "*TRG.bam" ! -name "*IGH.bam" ! -name "*IGK.bam" ! -name "*IGL.bam"`
			parallel -j${numCores} --eta  "samtools view -b -f4 {1} -h -o {1.}_unmapped.bam" ::: "${fileListM[@]}"
		fi
	fi
fi

##############################################################

if [[ "${toDo[*]}" == *extractReceptors* ]]; then
	echo "Extracting receptor regions..."
	echo `date +"%a %x %X"`
	regionTRA="${chr}14:20000000-24000000"
	regionTRB="${chr}7:140000000-145000000"
	regionTRD="${chr}14:20000000-25000000"
	regionTRG="${chr}7:36000000-41000000"
	regionIGH="${chr}14:103000000-107349000"
	regionIGK="${chr}2:86000000-92000000"
	regionIGL="${chr}22:21000000-25000000"

	fileListM=`find ${bamFolder} -type f -name "*.bam" ! -name "*unmapped*"  ! -name "*TRA.bam" ! -name "*TRB.bam" ! -name "*TRD.bam" ! -name "*TRG.bam" ! -name "*IGH.bam" ! -name "*IGK.bam" ! -name "*IGL.bam"`
	echo "Extracting TRA..."
	parallel -j${numCores} --eta  "samtools view -b -h {1} $regionTRA -o {1.}_TRA.bam" ::: "${fileListM[@]}"
	echo "Extracting TRB..."
	parallel -j${numCores} --eta  "samtools view -b -h {1} $regionTRB -o {1.}_TRB.bam" ::: "${fileListM[@]}"
	echo "Extracting TRD..."
	parallel -j${numCores} --eta  "samtools view -b -h {1} $regionTRD -o {1.}_TRD.bam" ::: "${fileListM[@]}"
	echo "Extracting TRG..."
	parallel -j${numCores} --eta  "samtools view -b -h {1} $regionTRG -o {1.}_TRG.bam" ::: "${fileListM[@]}"
	echo "Extracting IGH..."
	parallel -j${numCores} --eta  "samtools view -b -h {1} $regionIGH -o {1.}_IGH.bam" ::: "${fileListM[@]}"
	echo "Extracting IGK..."
	parallel -j${numCores} --eta  "samtools view -b -h {1} $regionIGK -o {1.}_IGK.bam" ::: "${fileListM[@]}"
	echo "Extracting IGL..."
	parallel -j${numCores} --eta  "samtools view -b -h {1} $regionIGL -o {1.}_IGL.bam" ::: "${fileListM[@]}"
fi

##############################################################
if [[ "${toDo[*]}" == *quickcheck* ]]; then
	echo "Running quickcheck..."
	echo `date +"%a %x %X"`
	numBams=`find $bamFolder -type f -name "*.bam" ! -name "*unmapped*"  ! -name "*TRA.bam" ! -name "*TRB.bam" ! -name "*TRD.bam" ! -name "*TRG.bam" ! -name "*IGH.bam" ! -name "*IGK.bam" ! -name "*IGL.bam" -printf '.' | wc -c`
	numBais=`find $bamFolder -type f -name "*.bai" -printf '.' | wc -c`
	
	echo "numBams = $numBams"
	echo "numBais = $numBais"

	if [ "$numBams" -ne "$numBais" ]; then
		echo "Warning: Some index files are missing. Indexing files..."
		echo `date +"%a %x %X"`

		fileListM=`find ${bamFolder} -type f -name "*.bam" ! -name "*unmapped*"  ! -name "*TRA.bam" ! -name "*TRB.bam" ! -name "*TRD.bam" ! -name "*TRG.bam" ! -name "*IGH.bam" ! -name "*IGK.bam" ! -name "*IGL.bam"`
		
		for bamFile in $fileListM
		do
			baseName=${bamFile%.bam}
			baiName="${baseName}.bai"

			if [ ! -e "$baiName" ]; then
				echo "Indexing ${bamFile}..."
				samtools index $bamFile
			fi
		done
	fi

	
	fileList=`find $bamFolder -type f -name "*.bam"`
	for bamFile in $fileList
	do
		samtools quickcheck -v $bamFile
	done
fi
##############################################################
if [[ "${toDo[*]}" == *runSearches* ]]; then
	echo "Running Searches..."
	echo `date +"%a %x %X"`

	fileListUM=`find ${bamFolder} -type f -name "*unmapped.bam"`
	fileListM=`find ${bamFolder} -type f -name "*.bam" ! -name "*unmapped*"  ! -name "*TRA.bam" ! -name "*TRB.bam" ! -name "*TRD.bam" ! -name "*TRG.bam" ! -name "*IGH.bam" ! -name "*IGK.bam" ! -name "*IGL.bam"`

	for IgTcr in "${IgTcrs[@]}"
	do
		echo "Processing $IgTcr..."
		echo `date +"%a %x %X"`
		if [[ "$IgTcr" == *UM* ]]; then
			parallel -j${numCores} --eta "sh $searchScript samtools {1} ${resultsDir}/${IgTcr}/{1/}.tsv ${IgTcr} ${chr}" ::: "${fileListUM[@]}"
		else
			parallel -j${numCores} --eta "sh $searchScript samtools {1} ${resultsDir}/${IgTcr}/{1/}.tsv ${IgTcr} ${chr}" ::: "${fileListM[@]}"
		fi
	done
fi
##############################################################

if [[ "${toDo[*]}" == *tsv2xlsx* ]]; then
	echo "Running tsv2xlsx..."
	echo `date +"%a %x %X"`

	for IgTcr in "${IgTcrs[@]}"
	do
		echo "Running tsv2xlsx ${IgTcr}..."
		echo `date +"%a %x %X"`

		mkdir -p "${resultsDir}_xlsx/$IgTcr"

		if [[ "$IgTcr" == *UM* ]]; then
			parallel -j${numCores} --eta "python $tsv2xlsx $resultsDir/${IgTcr}/{1/}.tsv ${resultsDir}_xlsx/${IgTcr}/{1/}.xlsx" ::: "${fileListUM[@]}"
		else
			parallel -j${numCores} --eta "python $tsv2xlsx $resultsDir/${IgTcr}/{1/}.tsv ${resultsDir}_xlsx/${IgTcr}/{1/}.xlsx" ::: "${fileListM[@]}"
		fi
	done
fi

##############################################################
if [[ "${toDo[*]}" == *imgtSearch* ]]; then
	echo "Running IMGT Search..."
	echo `date +"%a %x %X"`

	for IgTcr in "${IgTcrs[@]}"
	do
		echo "Running IMGT for ${IgTcr}..."
		echo `date +"%a %x %X"`
		$phpBin $imgtSearch $resultsDir/${IgTcr} $resultsDir/${samplesName}_${IgTcr}_vjMatchList.tsv ${IgTcr}
		if [[ "${toDo[*]}" == *tsv2xlsx* ]]; then
			python $tsv2xlsx $resultsDir/${samplesName}_${IgTcr}_vjMatchList.tsv ${resultsDir}_xlsx/${samplesName}_${IgTcr}_vjMatchList.xlsx
		fi
	done
fi

##############################################################

if [[ "${toDo[*]}" == *pigz* ]]; then
	tar -I pigz -cf ${resultsDir}.tar.gz ${resultsDir}
	if [[ "${toDo[*]}" == *tsv2xlsx* ]]; then
		tar -I pigz -cf ${resultsDir}_xlsx.tar.gz ${resultsDir}_xlsx
	fi
fi

##############################################################


if [[ "${toDo[*]}" == *emailEnd* ]]; then
	date +"%a %x %X" | mail -s "Script finished" $email
fi

##############################################################
echo "Finished."
echo `date +"%a %x %X"`