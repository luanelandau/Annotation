#!/bin/bash
#SBATCH --account=omergokc
#SBATCH --partition=general-compute
#SBATCH --qos=nih
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=800G
#SBATCH --job-name=GeMoMa_Opossum
#SBATCH --output=GeMoMa.Opossum1.out
#SBATCH --error=GeMoMa.Opossum1.err
#SBATCH --mail-user=luanejan@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --requeue
#Specifies that the job will be requeued after a node failure.
#The default is that the job will not be requeued.

#Gemoma (http://www.jstacs.de/index.php/GeMoMa#Frequently_asked_questions) is a software that uses another annotation to annotate the genome you're assemblyin. 
#The pool of tools of the Gemoma software goes far beyond that, in which Gemoma can use proteomes from multiple species or RNA-seq, so on and so forth. 
#I have not explored all the tools this software offers. Here, I am using the reference genome and its annotation for Opossum to annotate my own genome 
#using my own RNA-seq. This is on going work and I am still not sure of how the entire process goes.

#Defining variables
SPECIES="Opossum1" #Name of the species I want to annotate
TARGET="/projects/academic/omergokc/Luane/Opossum1/ragtag/ragtag_output/ragtag.scaffold.fasta" #My scaffold - not annotated
OUTDIR="/projects/academic/omergokc/Luane/Opossum1/Gemoma" #Directory I want to save my runs
GMM=/projects/academic/omergokc/Luane/conda/envs/annotation/share/gemoma-1.9-0 #path to the gemoma folder

#Details of the reference genome I am using to annotate mine
REF0="Monodelphis_domestica_GCF_027887165" #Just a string name for the reference genome used, can be whatever
REFGEN0="/projects/academic/omergokc/Luane/Opossum1/ncbi_dataset/data/GCF_027887165.1/GCF_027887165.1_mMonDom1.pri_genomic.fna" #path to the reference genome
REFASS0="/projects/academic/omergokc/Luane/Opossum1/ncbi_dataset/data/GCF_027887165.1/genomic.gff" #path to the annotation file

#My RNA-seq already mapped into my genome. This step has to be done earlier with hisat2. See previous scrips for this.
SAM_1="/projects/academic/omergokc/Luane/Opossum1/Opossum1_RNA/hisat2/Opos-01-BUC-Fem.sam"
SAM_2="/projects/academic/omergokc/Luane/Opossum1/Opossum1_RNA/hisat2/Opos-01-PAR-Fem-R.sam"
SAM_3="/projects/academic/omergokc/Luane/Opossum1/Opossum1_RNA/hisat2/Opos-01-SL-Fem-R.sam"
SAM_4="/projects/academic/omergokc/Luane/Opossum1/Opossum1_RNA/hisat2/Opos-01-SM-Fem-R.sam"

#Loading necessary modules
module load gcc/11.2.0 openmpi/4.1.1 gemoma/1.9 mmseqs2/13-45111

#Running gemoma. It requires you to use java, so make sure you have it. 
java -Xmx750g -jar $GMM/GeMoMa-1.9.jar CLI GeMoMaPipeline \
t=$TARGET o=true p=true \
s=own i=$REF0 a=$REFASS0 g=$REFGEN0 \
tblastn=false \
GeMoMa.c=0.4 \
GeMoMa.Score=ReAlign \
AnnotationFinalizer.r=NO \
Extractor.r=true \
Extractor.f=false \
r=MAPPED ERE.m=$SAM_1 ERE.m=$SAM_2 ERE.m=$SAM_3 ERE.m=$SAM_4 \
GAF.f="start=='M' and stop=='*' and (isNaN(score) or score/aa>='1.00')" \
outdir=$OUTDIR threads=32
