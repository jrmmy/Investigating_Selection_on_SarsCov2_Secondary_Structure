#!/bin/bash

#First argument = guidelines file
#Second argument = location where sequences are
#Third argument = location for output
#Fourth argument = baseline codeml ctl file

guidelinesFile=$1
inputFolder=$2
outputFolder=$3
basectl=$4

while read geneName minSize maxSize maxN
do
	echo $geneName
	geneOutputFolder=$outputFolder/${geneName}_output

	if [ ! -d "$geneOutputFolder" ]
	then
		mkdir $geneOutputFolder
	fi

	slurmFile="$geneOutputFolder/$geneName.slurm"
	echo "#!/bin/bash" > $slurmFile
	echo "module load TranslatorX" >> $slurmFile
	echo "module load seaview" >> $slurmFile
	echo "translatorx_vLocal.pl -i ${inputFolder}/${geneName}_sequences.fna -o $geneOutputFolder/$geneName" >> $slurmFile
	echo "module load IQ-Tree/2.0-rc1" >> $slurmFile
	echo "iqtree -s $geneOutputFolder/$geneName.nt_ali.fasta --prefix $geneOutputFolder/$geneName" >> $slurmFile
	ctlFile="$geneName.ctl"
	echo "sed 's|seqfile_to_replace|$geneOutputFolder/$geneName.nt_ali.fasta|' $basectl > $inputFolder/$ctlFile" >> $slurmFile
	echo "sed -i 's|treefile_to_replace|$geneOutputFolder/$geneName.tree|' $ctlFile" >> $slurmFile
	echo "sed -i 's|mlc_to_replace|$geneOutputFolder/$geneName.mlc|' $ctlFile" >> $slurmFile
#	echo "module load PAML" >> $slurmFile
#	echo "codeml $ctlFile" >> $slurmFile
	echo "eval \"\$(conda shell.bash hook)\"" >> $slurmFile 
	echo "conda activate /mnt/home/software/conda-envs/biopython/" >> $slurmFile
	echo "/mnt/home/jd2669/SarsCov2Project/Scripts/trim.py $geneOutputFolder/$geneName.nt_ali.fasta > $geneOutputFolder/$geneName.fna.trimmed" >> $slurmFile

	echo "conda activate /mnt/home/software/conda-envs/hyphy/" >> $slurmFile
	echo "hyphy slac $geneOutputFolder/$geneName.fna.trimmed --tree  $geneOutputFolder/${geneName}.treefile --samples 0 --pvalue 1 > ${geneOutputFolder}/${geneName}.slac" >> $slurmFile
	echo "grep \" p = \" ${geneOutputFolder}/${geneName}.slac | sed 's/ \+/,/g' | sed 's/|,//g' | sed 's/,|$//' | cut -d ',' -f 1,2,3,4,5,6,10 > ${geneOutputFolder}/${geneName}.slac.csv" >> $slurmFile
	#sbatch $slurmFile
	sleep 1s
done < $guidelinesFile




