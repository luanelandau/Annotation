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

#Making outdir
mkdir hisat2_refgenome
cd hisat2_refgenome
 
#Running hisat2
hisat2-build -f $ASSEMBLY Opossum_ref #Builds an index for the reference genome in hisat2 format

#maps each paired RNA-seq file to the corresponding genome using the index file generated in last step
hisat2 --threads 16 -x Opossum_ref -1 ${path}/${PAR}_RNA_R1.trimmed.fastq.gz -2 ${path}/${BUC}_RNA_R2.trimmed.fastq.gz -S ${BUC}_ref.sam
hisat2 --threads 16 -x Opossum_ref -1 ${path}/${PAR}_RNA_R1.trimmed.fastq.gz -2 ${path}/${PAR}_RNA_R2.trimmed.fastq.gz -S ${PAR}_ref.sam
hisat2 --threads 16 -x Opossum_ref -1 ${path}/${SL}_RNA_R1.trimmed.fastq.gz -2 ${path}/${SL}_RNA_R2.trimmed.fastq.gz -S ${SL}_ref.sam
hisat2 --threads 16 -x Opossum_ref -1 ${path}/${SM}_RNA_R1.trimmed.fastq.gz -2 ${path}/${SM}_RNA_R2.trimmed.fastq.gz -S ${SM}_ref.sam

#Organizing files to run featurecounts
module load samtools

# sorts the sam files & creates bam files
samtools sort hisat2/${PAR}_ref.sam  -O bam -o sortsam/${BUC}_ref.final.bam 
samtools sort hisat2/${PAR}_ref.sam  -O bam -o sortsam/${PAR}_ref.final.bam
samtools sort hisat2/${PAR}_ref.sam  -O bam -o sortsam/${SL}_ref.final.bam
samtools sort hisat2/${PAR}_ref.sam  -O bam -o sortsam/${SM}_ref.final.bam

#At this step, one can remove the sam files if wanted

# loads subread
module load subread

#counts genomic features such as exons, genes, promoters,
featureCounts -a $REFGFF -o readCounts/${BUC}.counts.txt -T 6 sortsam/${BUC}_ref.final.bam 
featureCounts -a $REFGFF -o readCounts/${PAR}.counts.txt -T 6 sortsam/${PAR}_ref.final.bam 
featureCounts -a $REFGFF -o readCounts/${SL}.counts.txt -T 6 sortsam/${SL}_ref.final.bam 
featureCounts -a $REFGFF -o readCounts/${SM}.counts.txt -T 6 sortsam/${SM}_ref.final.bam 
