#!/bin/bash
#SBATCH --account=omergokc
#SBATCH --partition=general-compute
#SBATCH --qos=nih
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --mem=400G
#SBATCH --ntasks-per-node=1
#SBATCH --job-name="cutadapt_fastqc_Opossum1"
#SBATCH --output=cutadapt_fastqc_Opossum1.out
#SBATCH --error=cutadapt_fastqc_Opossum1.err
#SBATCH --mail-user=luanejan@buffalo.edu
#SBATCH --mail-type=ALL
##SBATCH --requeue
#SBATCH --array=0-3 # job array index
#Specifies that the job will be requeued after a node failure.
#The default is that the job will not be requeued.

#This is a script to run cutadapt (https://github.com/marcelm/cutadapt) and fastqc (https://github.com/s-andrews/FastQC) for raw RNA-seq data. 
#Once you have your fastq files from illumina, you need to perform this processsing to use the RNA-seq data. 
#It will trimm the adaptors, delete reads below quality score of 20 (or whatever you set) and remove reads shorter than 15 bp (again, you can set to whatever). 
#Trim galore (https://github.com/FelixKrueger/TrimGalore) is another tool that uses both of these to perform the same things in one command line. 
#However, I like to do it step by step, so I decided to use these tools. 

#make a txt file with the prefixes of the fastq files you're processing. make sure the number of prefixes matched the arrays - 0 counts as 1 in this case. 
#Like:
#nano jobs_RNAseq.txt
#Opos-01-BUC-Fem
#Opos-01-PAR-Fem-R
#Opos-01-SL-Fem-R
#Opos-01-SM-Fem-R

#this will make it read each name from the txt file in array
names=($(cat jobs_RNAseq.txt))
path="/projects/academic/omergokc/Luane/Opossum1/Opossum1_RNA/" #path to where your RNA-seq is located

#make the directories for your outputs
mkdir trimmed fastqc

#load cutadapt and fastqc
module load gcccore/11.2.0 cutadapt/3.5 fastqc

#run cutadapt. This script is: trimming the adapters, filtering for reads with quality score more than 20, and minimum lenght of 15bp.
#These are the usual adaptors of illumina sequencing. Trim-galore can identify the adaptors automatically, which makes it user friendly in case you dont know what you're working with. 
#The --cores=0 option makes it use the cores defined in the script. But Im not sure it works well with the array option.
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -n 2 -o trimmed/${names[${SLURM_ARRAY_TASK_ID}]}_RNA_R1.trimmed.fastq.gz -p trimmed/${names[${SLURM_ARRAY_TASK_ID}]}_RNA_R2.trimmed.fastq.gz -q 20 --minimum-length 15 ${path}${names[${SLURM_ARRAY_TASK_ID}]}-RNA_R1_001.batch4.fastq.gz ${path}${names[${SLURM_ARRAY_TASK_ID}]}-RNA_R2_001.batch4.fastq.gz

#quality check on your RNA-seq files
fastqc trimmed/${names[${SLURM_ARRAY_TASK_ID}]}_RNA_R1.trimmed.fastq.gz -o fastqc/
fastqc trimmed/${names[${SLURM_ARRAY_TASK_ID}]}_RNA_R2.trimmed.fastq.gz -o fastqc/
