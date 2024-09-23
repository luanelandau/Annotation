#!/bin/bash
#SBATCH --account=omergokc
#SBATCH --partition=general-compute
#SBATCH --qos=nih
#SBATCH --time=2:00:00
#SBATCH --nodes=8
#SBATCH --mem=100G
#SBATCH --ntasks-per-node=1
#SBATCH --job-name="sra_download"
#SBATCH --output=SRA_dowload.out
#SBATCH --error=SRA_dowload.err
#SBATCH --mail-user=luanejan@buffalo.edu
#SBATCH --mail-type=ALL
##SBATCH --requeue
#SBATCH --array=0-26 # job array index
#Specifies that the job will be requeued after a node failure.
#The default is that the job will not be requeued.
#names=($(cat SRA_accession.txt))
names=($(cat files.txt))

eval "$(/projects/academic/omergokc/own_modules/2023.01/2023.01/software/Core/anaconda3/2023.03-1/bin/conda shell.bash hook)"
conda activate Assembly

module load gcc/11.2.0  openmpi/4.1.1 sra-toolkit

#prefetch ${names[${SLURM_ARRAY_TASK_ID}]}
fastq-dump --split-files ${names[${SLURM_ARRAY_TASK_ID}]}/*.sra 

#After, run this command on terminal to rename the files:
#while read -r srx srr name; do
#	mv ${srr}_1.fastq ${name}_R1.fastq
#done < SRA_accession_ids.txt
#Make sure to have a file SRA_accession_ids.txt that looks like this or something similar:

#SRX7566894 SRR10898068 Fetal_Parotid_F.475_S17
#SRX7566893 SRR10898067 Fetal_Submandibular_F.446_S6
#SRX7566892 SRR10898066 Fetal_Sublingual_F.445_S21
#SRX7566891 SRR10898065 Fetal_Submandibular_F.444_S23
#SRX7566890 SRR10898064 Fetal_Sublingual_F.431_S20
#SRX7566889 SRR10898063 Fetal_Submandibular_F.257_S5
#SRX7566888 SRR10898061 Fetal_Sublingual_F.255_S22
#SRX7566887 SRR10898060 Fetal_Sublingual_F.248_S18
#SRX7566886 SRR10898059 Fetal_Submandibular_F.22_S24
#SRX7566885 SRR10898058 Fetal_Parotid_F.21_S1
#SRX7566884 SRR10898057 Fetal_Submandibular_F.208_S4
#SRX7566883 SRR10898056 Fetal_Sublingual_F.202_S19
#SRX7566882 SRR10898055 Fetal_Parotid_F.165_S2
#SRX7566881 SRR10898054 Fetal_Sublingual_F.164_S3
#SRX7566880 SRR10898053 Adult_Sublingual_A.522_S16
#SRX7566879 SRR10898052 Adult_Submandibular_A.509_S14
#SRX7566878 SRR10898051 Adult_Submandibular_A.507_S13
#SRX7566877 SRR10898050 Adult_Submandibular_A.505_S26


