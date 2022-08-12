----------------------------------------
Lab Notebook - Demultiplexing Assignment
----------------------------------------
Logan Wallace 7/27/22
----------------------------------------------------------------------
Collaborators; Andreas M., Keenan R., Jack P.,Matt Esqueda, Peter Pham
----------------------------------------------------------------------

These are the files we are working with on talapas;
/projects/bgmp/shared/2017_sequencing/
1294_S1_L008_R1_001.fastq.gz
1294_S1_L008_R2_001.fastq.gz
1294_S1_L008_R3_001.fastq.gz
1294_S1_L008_R4_001.fastq.gz

Copied the bioinfo.py file to talapas
loganwallace$ scp bioinfo.py lwallac2@talapas-ln1.uoregon.edu:/projects/bgmp/lwallac2/bioinfo/Bi622/Demultiplex

Created testfile.fastq to peform practice runs with my python script before running the multi-hour script

- PART 1 - (find some bash commands below that were useful in answering each question)
1. i. R1 and R4 contain our biological sequences and R2/R3 contain our index reads. 
    $ ls -lah   initial data exploration, file size etc
    $ zcat 1294_S1_L008_R2_001.fastq.gz head -2 | tail -1 | wc
    $ zcat 1294_S1_L008_R2_001.fastq.gz | head -1000000 | tail -20      Just taking a look at 20 lines about 1,000,000 lines in
    $ wc -l 1294_S1_L008_R3_001.fastq.gz    How many lines are in this file?
    $ zcat 1294_S1_L008_R3_001.fastq.gz | grep -A1 "^@" | grep -v "^@" | grep -v "^--" | grep "N" | wc     Number of indexes containing 'N' 
    ii. The read lengths in each file are 101 for the biological reads and 8 for the index reads. 
    iii. The phred encoding is phred 33 for all four of the files. 
2. 'Generate a per base distribution of the quality scores for all four reads'.
    i. Need to generate 4 histograms for the quality scores. #1 Review python script from PS4 and run it on my testfile #2 Write an "sbatch" script to run a few jobs that will parse the enormous fastq files.
    ii. What is a good quality cutoff for index reads and biological read pairs to utilize for sample identification and downstream analysis. Commonly, Q30 is a standard for quality scores. This means an error rate of 1/1,000 bases called. 
    iii. 

Wrote both 'demultiplex_histogram.py' and 'run_demultiplex.sh'. 
The python script is going to run through each fastq file and return a histogram of the qscores.
The bash script creates a batch job for each of the fastq files
Outputs from this are named 'Read<x>.png'.

Note - The way I have written my test files; 
The first indexes should fail the quality filter and be placed into the unknown file. The second set should pass the filter and be placed in the match and the third have been written as hopped passing q filter, where both are in the index list but the second does not match the reverse complement of the other. 

- PART 2 -
Psuedocode
*See 'Psuedo' for my psuedocode.

- Part 3 -
*See 'Demultiplex.py' for my demultiplexing assignment.
Below, find the code for running Demultiplex.py assignment 
$ python Demultiplex.py -r1 read1_test.fastq -r2 read2_test.fastq -r3 read3_test.fastq -r4 read4_test.fastq

Wrote 'run_demux.sh'

Submitted sbatch for job ID 21973845 (fingers crossed); oof, everything wrote to the unknown file.

