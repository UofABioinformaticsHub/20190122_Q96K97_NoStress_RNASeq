#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=24:00:00
#SBATCH --mem=32GB
#SBATCH -o /fast/users/a1647910/20190122_Q96K97_NoStress_RNASeq/slurm/%x_%j.out
#SBATCH -e /fast/users/a1647910/20190122_Q96K97_NoStress_RNASeq/slurm/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lachlan.baer@adelaide.edu.au

# Define the number of cores
# This should match the selection above
CORES=16

# Define the project root
PROJROOT=/fast/users/a1647910/20190122_Q96K97_NoStress_RNASeq/

# Define any subdirectories
RAWFQ=${PROJROOT}/0_rawData/fastq
RAWQC=${PROJROOT}/0_rawData/FastQC

# Make sure all the directories exist, just by creating them
mkdir -p ${RAWQC}

# From here, we need to find all the files and then run fastqc

# Find the files
R1=$(ls ${RAWFQ}/*1.fq.gz)
echo -e "Found:\n${R1}"

# Now we need to step through this list and find the matching R2 file
for F1 in ${R1}
  do 
    echo "Currently analysing ${F1}"
    F2=${F1%1.fq.gz}2.fq.gz
    echo "The matching file should be ${F2}"
    # we shuold probably do a file.exists test on F2 here
    echo -e "
    fastqc \
      -t ${CORES} \
      --no-extract \
      -o ${RAWQC} \
      ${F1} ${F2} "

    exit

  done 
