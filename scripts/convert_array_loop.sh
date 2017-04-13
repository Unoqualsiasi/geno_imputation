#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --partition=hugemem,cigene
#SBATCH --output=slurm/convert_%A_%a.log  ##  when array job 
#SBATCH --mem-per-cpu=4000

set -o nounset   # Prevent unset variables from been used.
set -o errexit   # Exit if error occurs
nt=$SLURM_CPUS_ON_NODE ## sets $nt to -n antall cores requetsted

## Usage: sbatch -a 1-29 conform_25K_50K.sh arrayfile

arrayfile=$1
TASK=$SLURM_ARRAY_TASK_ID
infile=$(awk ' NR=='$TASK' { print $1 ; }' $arrayfile)
informat=$(awk ' NR=='$TASK' { print $2 ; }' $arrayfile)

# Usage: sh GenomeList_2_plink.sh [FinalReport file] [informat]

prefix=/mnt/users/tikn/for_folk/geno
sh $prefix/geno_imputation/scripts/Genome_2_plink.sh $infile $informat


