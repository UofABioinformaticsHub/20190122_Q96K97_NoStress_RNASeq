#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=4:00:00
#SBATCH --mem=32GB
#SBATCH -o /fast/users/a1647910/20190122_Q96K97_NoStress_RNASeq/slurm/%x_%j.out
#SBATCH -e /fast/users/a1647910/20190122_Q96K97_NoStress_RNASeq/slurm/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lachlan.baer@adelaide.edu.au

## Modules
module load picard/2.2.4-Java-1.8.0_71

## Directories
PROJROOT=/fast/users/a1647910/20190122_Q96K97_NoStress_RNASeq
ALIGNDATA=${PROJROOT}/2_alignedData
MARKDUPLICATES=${PROJROOT}/5_markDuplicates
PICARD=/apps/software/picard/2.2.4-Java-1.8.0_71/picard.jar

## Mark Duplicates
for i in ${ALIGNDATA}/bam/*.bam
  do

    echo -e "Currently working on ${i}"

    ## Create output name
    OUT=${MARKDUPLICATES}/$(basename $i)
    echo -e "The output file will be ${OUT}"
    ## Create log name
    LOG=${MARKDUPLICATES}/$(basename ${i%Aligned.sortedByCoord.out.bam})_metrics.txt
    echo -e "The log file will be ${LOG}"

    ## Run Picard
    java -jar $PICARD MarkDuplicates \
      I=${i} \
      O=${OUT} \
      M=${LOG} \
      VALIDATION_STRINGENCY=SILENT \
      CREATE_INDEX=TRUE

  done
