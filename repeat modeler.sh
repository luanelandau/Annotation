#!/bin/bash
#SBATCH --account=omergokc
#SBATCH --partition=general-compute
#SBATCH --qos=nih
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G
#SBATCH --job-name="RModeler_Opossum1"
#SBATCH --output=RModeler_Opossum1.out
#SBATCH --error=RModeler_Opossum1.err
#SBATCH --mail-user=luanejan@buffalo.edu
#SBATCH --mail-type=ALL

#This script if a first step into looking at repeated sequences on the genome, and further masking it. 
#It uses repeatmodeler to do this (https://github.com/Dfam-consortium/RepeatModeler/).

ASSEMBLY="ragtag.scaffold.fasta" #the assembly has to be on the same folder

module load gcc/11.2.0 openmpi/4.1.1 repeatmodeler/2.0.4.KRAB #loading the modules necessary

BuildDatabase -name Opossum1 $ASSEMBLY
RepeatModeler -database Opossum1 -threads 16 -LTRStruct

echo "Complete"
