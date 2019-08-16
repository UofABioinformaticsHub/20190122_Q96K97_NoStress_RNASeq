#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=1:00:00
#SBATCH --mem=32GB
#SBATCH -o /fast/users/a1647910/20190122_Q96K97_NoStress_RNASeq/slurm/%x_%j.out
#SBATCH -e /fast/users/a1647910/20190122_Q96K97_NoStress_RNASeq/slurm/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lachlan.baer@adelaide.edu.au

## Cores
CORES=16

########################################
## Full run of the data from Novogene ##
########################################

## Modules
module load FastQC/0.11.7
module load STAR/2.5.3a-foss-2016b
module load SAMtools/1.3.1-GCC-5.3.0-binutils-2.25
module load AdapterRemoval/2.2.1-foss-2016b
module load GCC/5.4.0-2.26
module load Subread

## Genomic Data Files
REFS=/fast/users/a1647910/20190122_Q96K97_NoStress_RNASeq/gtf_temp/
GTF=${REFS}Danio_rerio.GRCz11.94.chr.gtf

## # The STAR genome index was generated manually on clarence after mounting /data/biorefs
## # STAR 2.6.0a was used for generation of the index, whilst STAR 2.5.3 will be
## # used for alignment
## STAR \
##   --runThreadN 20 \
##   --runMode genomeGenerate \
##   --genomeDir ${REFS}/star/ \
##   --genomeFastaFiles ${REFS}/Danio_rerio.GRCz11.dna.primary_assembly.fa \
##   --sjdbGTFfile ${REFS}/Danio_rerio.GRCz11.94.chr.gtf \
##   --sjdbOverhang 74

## Directories
PROJROOT=/fast/users/a1647910/20190122_Q96K97_NoStress_RNASeq

## Directories for Initial FastQC
RAWDATA=${PROJROOT}/0_rawData
mkdir -p ${RAWDATA}/FastQC

## Setup for Trimmed data
TRIMDATA=${PROJROOT}/1_trimmedData
mkdir -p ${TRIMDATA}/fastq
mkdir -p ${TRIMDATA}/FastQC
mkdir -p ${TRIMDATA}/log

## Setup for genome alignment
ALIGNDATA=${PROJROOT}/2_alignedData
mkdir -p ${ALIGNDATA}/log
mkdir -p ${ALIGNDATA}/bam
mkdir -p ${ALIGNDATA}/FastQC
mkdir -p ${ALIGNDATA}/featureCounts


##--------------------------------------------------------------------------------------------##
## FastQC on the raw data
##--------------------------------------------------------------------------------------------##

#fastqc -t ${CORES} -o ${RAWDATA}/FastQC --noextract ${RAWDATA}/fastq/*.fq.gz

##--------------------------------------------------------------------------------------------##
## Trimming the Merged data
##--------------------------------------------------------------------------------------------##

#for R1 in ${RAWDATA}/fastq/*1.fq.gz
#  do
#
#    echo -e "Currently working on ${R1}"
#
#    # Now create the output filenames
#    out1=${TRIMDATA}/fastq/$(basename $R1)
#    BNAME=${TRIMDATA}/fastq/$(basename ${R1%_1.fq.gz})
#    R2=${R1%_1.fq.gz}_2.fq.gz
#    out2=${TRIMDATA}/fastq/$(basename $R2)
#    echo -e "The R2 file should be ${R2}"
#    echo -e "Output file 1 will be ${out1}"
#    echo -e "Output file 2 will be ${out2}"
#    echo -e "Trimming:\t${BNAME}"
#
#    #Trim
#    AdapterRemoval \
#      --gzip \
#      --trimns \
#      --trimqualities \
#      --minquality 30 \
#      --minlength 35 \
#      --threads ${CORES} \
#      --basename ${BNAME} \
#      --output1 ${out1} \
#      --output2 ${out2}\
#      --file1 ${R1} \
#      --file2 ${R2}
#
#  done

## Move the log files into their own folder
#mv ${TRIMDATA}/fastq/*settings ${TRIMDATA}/log

## Run FastQC
#fastqc -t ${CORES} -o ${TRIMDATA}/FastQC --noextract ${TRIMDATA}/fastq/*.fq.gz


##--------------------------------------------------------------------------------------------##
## Aligning trimmed data to the genome
##--------------------------------------------------------------------------------------------##

## Aligning, filtering and sorting
#for R1 in ${TRIMDATA}/fastq/*1.fq.gz
#  do
#
#  BNAME=$(basename ${R1%_1.fq.gz})
#  R2=${R1%_1.fq.gz}_2.fq.gz
#  echo -e "STAR will align:\t${R1}"
#  echo -e "STAR will also align:\t${R2}"
#
#    STAR \
#        --runThreadN ${CORES} \
#        --genomeDir ${REFS}/star \
#        --readFilesIn ${R1} ${R2} \
#        --readFilesCommand gunzip -c \
#        --outFileNamePrefix ${ALIGNDATA}/bam/${BNAME} \
#        --outSAMtype BAM SortedByCoordinate
#
#  done

## Move the log files into their own folder
#mv ${ALIGNDATA}/bam/*out ${ALIGNDATA}/log
#mv ${ALIGNDATA}/bam/*tab ${ALIGNDATA}/log
#
# Fastqc and indexing
#for BAM in ${ALIGNDATA}/bam/*.bam
# do
#   fastqc -t ${CORES} -f bam_mapped -o ${ALIGNDATA}/FastQC --noextract ${BAM}
#   samtools index ${BAM}
# done

##--------------------------------------------------------------------------------------------##
## featureCounts
##--------------------------------------------------------------------------------------------##

## Feature Counts - obtaining all sorted bam files
sampleList=`find ${ALIGNDATA}/bam -name "*out.bam" | tr '\n' ' '`

## Running featureCounts on the sorted bam files
featureCounts -Q 10 \
  -s 2 \
  -T ${CORES} \
  -p \
  -a ${GTF} \
  -o ${ALIGNDATA}/featureCounts/counts.out ${sampleList}

## Storing the output in a single file
cut -f1,7- ${ALIGNDATA}/featureCounts/counts.out | \
sed 1d > ${ALIGNDATA}/featureCounts/genes.out
