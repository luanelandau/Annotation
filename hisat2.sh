#!/bin/bash
#SBATCH --account=omergokc
#SBATCH --partition=general-compute
#SBATCH --qos=nih
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=400G
#SBATCH --job-name="hisat_Opossum1"
#SBATCH --output=hisat_Opossum1.out
#SBATCH --error=hisat_Opossum1.err
#SBATCH --mail-user=luanejan@buffalo.edu
#SBATCH --mail-type=ALL
##SBATCH --requeue
#Specifies that the job will be requeued after a node failure.
#The default is that the job will not be requeued.

#HISAT2 is the main software used for mapping RNA-seq files to assemblies. This step is necessary for annotating RNA-seq files or for de novo assembly. 

#These are the prefixes of my RNA-seq files.
BUC="Opos-01-BUC-Fem"
PAR="Opos-01-PAR-Fem-R"
SL="Opos-01-SL-Fem-R"
SM="Opos-01-SM-Fem-R"

path="/projects/academic/omergokc/Luane/Opossum1/Opossum1_RNA/trimmed" #path of where the RNA-seq files are

ASSEMBLY="/projects/academic/omergokc/Luane/Opossum1/ragtag/ragtag_output_hap1/ragtag.scaffold.fasta"

module load gcc/11.2.0 openmpi/4.1.1 hisat2/2.2.1

#Make and name the directory you want the outputs to be
#mkdir hisat2
#Open the directory
cd hisat2

#Have to build the assembly index using hisat
hisat2-build -f $ASSEMBLY Opossum1

#This step will map in pairs the RNA-seq to the genome and generate .sam files. 
hisat2 --threads 16 -x Opossum1 -1 ${path}/${PAR}_RNA_R1.trimmed.fastq.gz -2 ${path}/${PAR}_RNA_R2.trimmed.fastq.gz -S ${PAR}.sam

hisat2 --threads 16 -x Opossum1 -1 ${path}/${SL}_RNA_R1.trimmed.fastq.gz -2 ${path}/${SL}_RNA_R2.trimmed.fastq.gz -S ${SL}.sam

hisat2 --threads 16 -x Opossum1 -1 ${path}/${SM}_RNA_R1.trimmed.fastq.gz -2 ${path}/${SM}_RNA_R2.trimmed.fastq.gz -S ${SM}.sam
