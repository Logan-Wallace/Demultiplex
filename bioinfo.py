# Author: Logan Wallace lwallac2@uoregon.edu

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.'''

__version__ = "0.3"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning



DNA_bases = "ACTG"
RNA_bases = "ACUG"



def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    phred: int = ord(letter) - 33
    return(phred)



def qual_score():
    ''''''
    pass



#Ask Leslie about flags;
def validate_base_seq(sequence, RNA_flag):
    '''This function takes a string as input and will return True if the string is DNA/RNA and False otherwise'''
    sequence = sequence.upper()
    return len(sequence) == sequence.count("A") + sequence.count("C") + sequence.count("G") +  sequence.count("U" if RNAflag else "T")



def gc_content():
    '''This function will take a string as an argument and return the percent of G's and C's as a whole number'''
    assert validate_DNA_seq(DNA), "Calling validate_base_seq above to ensure that we have a sequence of base pairs"
    DNA = DNA.upper()
    G = DNA.count("G")
    C = DNA.count("C")
    return ((G+C)/len(DNA))*100



def oneline_fasta(read_file: str):
    '''This function will take as input a fasta file and write out to a file called "oneline_fasta"'''
    '''Note that this function operates by using a string slicing method (see 'if line[] loop') that is specific to each file!!!'''
    #init some variables
    read_file: str = ""
    gene_ID: str = ""
    peptide: str = ""
    protein_ID: str = ""
    gene_name: str = ""


    #init a dict to hold our genes and sequences from the fasta
    peptide_dict:dict = {}

    #open our fasta file to read
    with open(read_file, 'r') as fasta:
        #loop through our fasta file and grab the gene_stable_ID
        for line in fasta:
            #if we are at a header line grab the gene ID
            if line[0] == '>':
                gene = line.split(" ")
                gene_ID = gene[3]
                gene_ID = gene_ID[5:20]
                protein_ID = gene[0]
                protein_ID = protein_ID[1:17]
                try: 
                    gene_name = gene[7]
                    gene_name = gene_name[12:]
                    gene_name = gene_name.strip("\n")
                except:
                    gene_name = "NA"
            #if not a header line grab the peptide sequence as 'peptide'
            else:
                peptide += line
                peptide = peptide.strip("\n")
                #if the gene is in the dictionary
                if gene_ID in peptide_dict:
                    #and the length of the peptide sequence is greater, add to dictionary
                    if len(peptide_dict[gene_ID][2]) > len(peptide):
                        peptide_dict[gene_ID] = (protein_ID, gene_name, peptide)
                #if the gene isn't in the dictionary, 
                #update with the gene and peptide sequence
                else:
                    peptide_dict[gene_ID] = (protein_ID, gene_name, peptide)
                    peptide = ""

    #open the file to read and one to write to new
    with open("oneline_fasta", 'w') as new:
        for i in peptide_dict:
            new.write(">" + peptide_dict[i][0] + " " + i + " " + peptide_dict[i][1] + "\n" + peptide_dict[i][2] + "\n")



def init_list(lst: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    for n in range(101):
        lst.append(value)
    return lst



def populate_list(file: str):
    #This function will return an array and a counter that will 
    #store quality score sums for each position of a read in a fastq file
    #create an empty list with init list
    q_list: list = []
    q_list = init_list(q_list)
    lcount: int = 0
    # #import the gzip module so we can read zipped files without unzipping them previously
    # import gzip
    #open our fastq file as 'fh'   
    with open(file, "r") as fh:
    #for all the lines in fh 
        for line in fh:
        #increment a line counter for each line read
            lcount += 1
        #if it is a phred score line
            if lcount % 4 == 0:
            #for each character in the line
                for i in range(len(line.strip())):
                #add the value of the phred score to that position in q_list
                    q_list[i] += convert_phred(line[i])
    return (q_list, lcount)



def knorm(file: str):
    '''this function will take in a fastq file as input and others (see arguments below) and return
    a normalized kmer fastq file.'''
    #init my global variables
    import argparse
    parser = argparse.ArgumentParser(description = "Arguments for kspec")
    parser.add_argument("-c", "--coverage_limit", help = "", type = int, default = 10)
    parser.add_argument("-f", "--filename", help = "", type = str, default = "lane1_NoIndex_L001_R1_003.fastq")
    parser.add_argument("-k", "--kmer_length", help = "", type = int, default = 15)
    parser.add_argument("-o", "--output_file", help = "", type = str, default = "knorm_l1_10x")
    args = parser.parse_args()

    #Below I am making it so that I can reference my variables within my file and 
    #avoid having to call 'args.' each time.
    c = args.coverage_limit
    f = args.filename
    k = args.kmer_length
    o = args.output_file

    #init an empty dict 'kmer_occur' to store 
    #keys=kmers, values=occurrences
    kmer_occur = {}

    #open a file to output reads to keep 
    with open(o, 'w') as out:
        #init a counter to keep track of lines in the fastq file
        counter = 0
        #open our fastq file 
        with open(f, 'r') as fastq:
            #loop through our fastq and record
            #each line to a variable
            for line in fastq:
                line = line.strip("\n")
                counter +=1
                if counter%4==1:
                    a = (line)
                elif counter%4==2:
                    b = (line)
                elif counter%4==3:
                    c1 = (line)            
                elif counter%4==0:
                    d = (line)  
                #init a list to store kmer coverage
                    kmer_coverage=[]
                    #for the length of the read
                    for x in range(len(b) - k + 1):
                        #slice the sequence by the length of a kmer
                        e =(b[x:x+k])
                        #if mer not in kmer_occur add to the dictionary with value 0
                        if e in kmer_coverage:
                            kmer_occur[e]+=1
                        #if mer in kmer_occur increment the value by one 
                        else: 
                            kmer_occur[e]=1
                    for y in range(len(b) - k+1):
                        e = (b[y:y+k])
                        kmer_coverage.append(kmer_occur[e])
                    kmer_coverage.sort()

                    position = int(len(kmer_occur)//2)
                    kmer_coverage_median = kmer_coverage[position]
                    if kmer_coverage_median < c:
                        final_file=(a + "\n" + b + "\n" + c1 + "\n" + d + "\n")
                        out.write(final_file)
                    kmer_occur={}
    return{kmer_occur}



if __name__ == "__main__":
    '''below are unit tests for some of the above functions'''
    #gc_content
    assert gc_content("GCGCGC") == 100
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATGCAT") == 50
    print("correctly calculated GC content")
    #validate_base_seq
    assert validate_DNA_seq("aaaaa") == True, "DNA string not recognized"
    print("Correctly identified a DNA string")
    assert validate_DNA_seq("Hi there!") == False, "Non-DNA identified as DNA"
    print("Correctly determined non-DNA")