#!/bin/bash
#SBATCH --account=bgmp                 
#SBATCH --partition=bgmp               
#SBATCH --job-name=Demultiplex                           
#SBATCH --nodes=1                      
#SBATCH --cpus-per-task=8


#python demultiplex.py -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz" -o "Read1" -r 101
#python demultiplex.py -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz" -o "Read2" -r 8
#python demultiplex.py -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz" -o "Read3" -r 8
python demultiplex.py -f "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz" -o "Read4" -r 101
