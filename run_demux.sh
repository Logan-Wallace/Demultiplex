#!/bin/bash
#SBATCH --account=bgmp                 
#SBATCH --partition=bgmp               
#SBATCH --job-name=Demultiplex         
#SBATCH --nodes=1                   
#SBATCH --cpus-per-task=1      

conda activate base
python Demultiplex.py 