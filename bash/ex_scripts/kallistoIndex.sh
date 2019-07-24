#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=01:00:00
#SBATCH --mem=4GB
#SBATCH -o /fast/users/a1647910/20190122_Q96K97_NoStress_RNASeq/slurm/%x_%j.out
#SBATCH -e /fast/users/a1647910/20190122_Q96K97_NoStress_RNASeq/slurm/%x_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lachlan.baer@adelaide.edu.au

#Load modules
module load kallisto/0.43.1-foss-2016b
module load SAMtools/1.3.1-foss-2016b

#Run kallisto
kallisto index \
  -i /fast/users/a1647910/20190122_Q96K97_NoStress_RNASeq/refs/Danio_rerio.GRCz11.cdna.inc.201psen1 \
  /fast/users/a1647910/20190122_Q96K97_NoStress_RNASeq/refs/Danio_rerio.GRCz11.cdna.inc.201psen1.fa
