----------------------------------------
Lab Notebook - Demultiplexing Assignment
----------------------------------------
Logan Wallace 7/27/22
----------------------------------------------
Collaborators; Andreas M., Keenan R., Jack P.,
----------------------------------------------

These are the files we are working with on talapas;
/projects/bgmp/shared/2017_sequencing/
1294_S1_L008_R1_001.fastq.gz
1294_S1_L008_R2_001.fastq.gz
1294_S1_L008_R3_001.fastq.gz
1294_S1_L008_R4_001.fastq.gz

Copied the bioinfo.py file to talapas
loganwallace$ scp bioinfo.py lwallac2@talapas-ln1.uoregon.edu:/projects/bgmp/lwallac2/bioinfo/Bi622/Demultiplex

Created testfile.fastq to peform practice runs with my python script before running the multi-hour script

- PART 1 -
1. i. R1 and R4 contain our biological sequences and R2/R3 contain our index reads. 
    ii. The read lengths in each file are 101 for the biological reads and 8 for the index reads. 
    iii. The phred encoding is phred 33 for all four of the files. 
2. 'Generate a per base distribution of the quality scores for all four reads'.
    i. Need to generate 4 histograms for the quality scores. #1 Review python script from PS4 and run it on my testfile #2 Write an "sbatch" script to run a few jobs that will parse the enormous fastq files.
    ii. What is a good quality cutoff for index reads and biological read pairs to utilize for sample identification and downstream analysis. Commonly, Q30 is a standard for quality scores. This means an error rate of 1/1,000 bases called. 
    iii. 

Wrote both 'demultiplex.py' and 'run_demultiplex.sh'. 
The python script is going to run through each fastq file and return a histogram of the qscores.
The bash script creates a batch job for each of the fastq files
Outputs from this are named 'Readx.png'.

- PART 2 -
Psuedocode
