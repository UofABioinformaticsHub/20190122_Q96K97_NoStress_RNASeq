
#!/bin/bash

## Directories
PROJROOT=/fast/users/a1647910/20190122_Q96K97_NoStress_RNASeq
TRIMDATA=${PROJROOT}/1_trimmedData

## Setup for kallisto output
mkdir -p ${PROJROOT}/3_kallisto

## Fire off the alignments
FQ=$(ls ${TRIMDATA}/fastq/*1.fq.gz)
echo -e "Found:\n${FQ}"

for F1 in ${FQ}
	do
	sbatch ${PROJROOT}/bash/kallisto.sh ${F1}
done
