# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | Biological Read | 101 | 33 |
| 1294_S1_L008_R2_001.fastq.gz | Index Read | 8 | 33 |
| 1294_S1_L008_R3_001.fastq.gz | Index Read | 8 | 33 |
| 1294_S1_L008_R4_001.fastq.gz | Biological Read | 101 | 33 |

2. Per-base NT distribution
    i. Use markdown to insert your 4 histograms here.
    
    #Read 2
    ![Read2](https://user-images.githubusercontent.com/107602241/181381885-34874444-e5a8-485c-ac10-a0623de47bd3.png)
    
    "Read 3"
    ![Read3](https://user-images.githubusercontent.com/107602241/181381905-73152b12-b336-4120-8965-90df199fb039.png)
    
    ii. Q30 is a standard cutoff score when sequencing and is also suggested to be a good minimum by Illumna 
    themselves. 
    It would mean an error rate in base calls of 1/1,000. Between coverage and the fact that our biological reads are 
    only 101bp long I would suggest
    this as a good standard for the biological reads. In the end we are simply trying to align the reads to a genome. 
    For our barcodes, the reads are more than 10x shorter so errors should be less common within a 
    barcode read than a biological read. However, it is very important that our barcodes are correct as this is how we 
    'demultiplex our samples
    after sequencing. I would suggest a cutoff of 33 for our barcodes (for each base) as this is quite stringent and 
    should help prevent barcode hopping. This is also based upon our results. The mean q-score for the first few reads 
    is just under this cutoff so we risk ditching a lot of data right away.
    iii. zcat 1294_S1_L008_R2_001.fastq.gz | grep -A1 "^@" | grep -v "^@" | grep -v "^--" | grep "N" | wc
         R2 - 3,976,613
         R3 - 3,328,051
    
## Part 2
1. Define the problem:
2. Describe output:
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ):
4. Pseudocode:
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
