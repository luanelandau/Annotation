#!/bin/bash
#SBATCH --account=omergokc
#SBATCH --partition=general-compute
#SBATCH --qos=nih
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=150G
#SBATCH --job-name="hisat_Opossum1"
#SBATCH --output=hisat_Opossum1.out
#SBATCH --error=hisat_Opossum1.err
#SBATCH --mail-user=luanejan@buffalo.edu
#SBATCH --mail-type=ALL


ASSEMBLY="/projects/academic/omergokc/Luane/Opossum1/ragtag/ragtag_output_hap1/ragtag.scaffold.fasta"

#The following are my RNA-seq files. I am just going to unzip them all because I need fastq files for this analysis. So next two steps are about that.
#RNA_SEQ_PAR_R1="/projects/academic/omergokc/Luane/Opossum1/Opossum1_genome_annotation/Opossum1_RNA/Opos-01-PAR-Fem-R-RNA_R1_001.batch4.fastq.gz"
#RNA_SEQ_PAR_R2="/projects/academic/omergokc/Luane/Opossum1/Opossum1_genome_annotation/Opossum1_RNA/Opos-01-PAR-Fem-R-RNA_R2_001.batch4.fastq.gz"
#RNA_SEQ_SL_R1="/projects/academic/omergokc/Luane/Opossum1/Opossum1_genome_annotation/Opossum1_RNA/Opos-01-SL-Fem-R-RNA_R1_001.batch4.fastq.gz"
#RNA_SEQ_SL_R2="/projects/academic/omergokc/Luane/Opossum1/Opossum1_genome_annotation/Opossum1_RNA/Opos-01-SL-Fem-R-RNA_R2_001.batch4.fastq.gz"
#RNA_SEQ_SM_R1="/projects/academic/omergokc/Luane/Opossum1/Opossum1_genome_annotation/Opossum1_RNA/Opos-01-SM-Fem-R-RNA_R1_001.batch4.fastq.gz"
#RNA_SEQ_SM_R2="/projects/academic/omergokc/Luane/Opossum1/Opossum1_genome_annotation/Opossum1_RNA/Opos-01-SM-Fem-R-RNA_R2_001.batch4.fastq.gz"
#RNA_SEQ_BUC_R1="/projects/academic/omergokc/Luane/Opossum1/Opossum1_genome_annotation/Opossum1_RNA/Opos-01-BUC-Fem-RNA_R1_001.batch4.fastq.gz"
#RNA_SEQ_BUC_R2="/projects/academic/omergokc/Luane/Opossum1/Opossum1_genome_annotation/Opossum1_RNA/Opos-01-BUC-Fem-RNA_R2_001.batch4.fastq.gz"

#gunzip -c $RNA_SEQ_PAR_R1 > Opos-01-PAR_R1.fastq
#gunzip -c $RNA_SEQ_PAR_R2 > Opos-01_PAR_R2.fastq
#gunzip -c $RNA_SEQ_SL_R1 >  Opos-01_SL_R1.fastq
#gunzip -c $RNA_SEQ_SL_R2 >  Opos-01_SL_R2.fastq
#gunzip -c $RNA_SEQ_SM_R1 >  Opos-01_SM_R1.fastq
#gunzip -c $RNA_SEQ_SM_R2 >  Opos-01_SM_R2.fastq
#gunzip -c $RNA_SEQ_BUC_R1 > Opos-01_BUC_R1.fastq
#gunzip -c $RNA_SEQ_BUC_R2 > Opos-01_BUC_R2.fastq

#I had a ridiculous problem doing this because I left a space after the "=". Now it should be all working

PAR_R1="Opos-01-PAR_R1.fastq"
PAR_R2="Opos-01_PAR_R2.fastq"
SL_R1="Opos-01_SL_R1.fastq"
SL_R2="Opos-01_SL_R2.fastq"
SM_R1="Opos-01_SM_R1.fastq"
SM_R2="Opos-01_SM_R2.fastq"
BUC_R1="Opos-01_BUC_R1.fastq"
BUC_R2="Opos-01_BUC_R2.fastq"

module load gcc/11.2.0 openmpi/4.1.1 hisat2/2.2.1

hisat2-build -f $ASSEMBLY Opossum1
hisat2 -x Opossum1 -1 $PAR_R1,$SL_R1,$SM_R1,$BUC_R1 -2 $PAR_R2,$SL_R2,$SM_R2,$BUC_R2 -S Opossum1_output.sam
