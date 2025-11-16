#!/bin/bash

#SBATCH --tasks 1
#SBATCH --cpus-per-task 1

# eval "$(conda shell.bash hook)"
conda activate /mnt/home/software/conda-envs/biopython/

# For other jobs
python pairwiseAlignment.py -i /mnt/scratch/goutlab/Virus_dNdS/jd2669/Batch10k/10kWholeGenome -o /mnt/scratch/goutlab/Virus_dNdS/jd2669/Batch10k/aligned10k -g /mnt/scratch/goutlab/Virus_dNdS/Data/geneGuide -f /mnt/scratch/goutlab/Virus_dNdS/jd2669/Batch10k/failed10k
# iqtree -s /mnt/scratch/goutlab/Virus_dNdS/jd2669/Batch5/batch5Chop_ValidWholeGenomes --prefix batch5tree


# For hyphy

#First argument = full path & prefix to the sequences
#Second argument = location for output

# directory=$1
# outputFolder=$2

# for file in "$directory"*; do
#     geneName=`basename $file`
#     echo "$file -> $geneName"
# 	  slurmFile="$outputFolder/$geneName.slurm"
#     #echo $slurmFile
#     echo "#!/bin/bash" > $slurmFile
#     echo "eval \"\$(conda shell.bash hook)\"" >> $slurmFile 
#     echo "conda activate /mnt/home/software/conda-envs/hyphy/" >> $slurmFile
#     echo "hyphy slac --alignment $file --tree /mnt/scratch/goutlab/Virus_dNdS/jd2669/Batch5/batch5tree.treefile --samples 0 --pvalue 1 > $outputFolder/$geneName-hyphy.log 2> $outputFolder/$geneName-hyphy.err" >> $slurmFile
#     sbatch $slurmFile
#     sleep 1s
     
# done 
