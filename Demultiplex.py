#! usr/env/bin python

#Logan Wallace 8/3/22

'''The purpose of this program is to parse through four different fastq files and demultiplex them. 
It will return a total of 52 files, 48 paired reads, 2 unknown and 2 hopped files. 
The program will also give metrics on the number of matched pairs between each combination of indices.'''

#import necessary modules
import argparse
import bioinfo
import itertools
import matplotlib as plt 
import gzip as gz

parser = argparse.ArgumentParser(description = "Arguments for demultiplex.py")
parser.add_argument("-r1", "--read_1", help = "The name of the fastq file for (biological) read1", type =str, default = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz")
parser.add_argument("-r2", "--read_2", help = "The name of the fastq file for (index) read2", type =str, default = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz")
parser.add_argument("-r3", "--read_3", help = "The name of the fastq file for (index) read3", type =str, default = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz")
parser.add_argument("-r4", "--read_4", help = "The name of the fastq file for (biological) read4", type =str, default = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz")
parser.add_argument("-i", "--index", help = "The name of the file containing all the indexes", type =str, default = "/projects/bgmp/shared/2017_sequencing/indexes.txt")

args = parser.parse_args()

#Variable initialization
index_list: list = []
index_dict: dict = {}
i = args.index
minimum_index_score: int = 25
minimum_biological_read_average: int = 0
index_pair_instances: dict = {}
r1 = args.read_1
r2 = args.read_2
r3 = args.read_3
r4 = args.read_4
matched_count: int = 0
unknown_count: int = 0
hopped_count: int = 0

#Functions
def reverse_complement(index: str):
    '''This function will return the reverse complement to the input string'''
    rcomp = ""
    for base in index:
        if base == "A":
            base = "T"
        elif base == "T":
            base = "A"
        elif base == "C":
            base = "G"
        else:
            base = "C"
        rcomp = base + rcomp
    #test: ATCG == CGAT
    return(rcomp)

def index_quality_check(index, index_qscores):
    '''This function will take in a string of index read and a string of phred 33 quality scores from index reads as input and determine if the string 
    meets the minimum quality scores or contains an "N" and then return this as either False (contains "N"), True (pass) or False (fail)'''
    #Does the index contain any N's?
    if "N" in index:
        return False 
    #Does the index fail the q-check at any point?
    for value in index_qscores:
        score = bioinfo.convert_phred(value)
        if score < minimum_index_score:
            #print(score, value)
            return False
    return True

def biological_quality_check(biological_qscores):
    '''This function will take in a string of quality score characters for a biological read and average them to make sure they meet a minimum of 30
    returning True if the minimum is met and False otherwise'''
    sum: int = 0
    for phred in biological_qscores:
        q = bioinfo.convert_phred(phred)
        sum += q
    average = sum / 101
    if average >= minimum_biological_read_average:
        return True
    else: 
        return False

def matched_hopped(index2, index3):
    '''This function will take in the index read sequences from reads 2 and 3 and tell us whether they are matched or hopped'''
    index3 = reverse_complement(index3)
    if index2 == index3:
        return True
    else:
        return False  

#Create a list of the indexes and a list of reverse complements of those indexes
with open(i) as indexes:
    count = 0
    for index in indexes:
        if count > 0:
            index = index.strip("\n").split("\t")[4]
            index_list.append(index)
        count += 1

#Initialize a dictionary with using index_list to later count on 
for index in index_list:
    index_dict[index] = 0

#Create a dictionary to hold instances of index pairs
#Remember to add in the self matches
index_combinations = list(itertools.permutations(index_list, r=2))
for combo in index_combinations:
    index_pair_instances.update({combo:0})
for index in index_list:
    index_pair_instances.update({(index, index):0})
#print(index_list)
#print(len(index_pair_instances))

r1_filenames_dict = {}
r2_filenames_dict = {}

#Generate dictionaries for r1/r2 filenames
for i in index_list:
    r1_filename = open("read1_" + i + ".fastq", "w") 
    r2_filename = open("read2_" + i + ".fastq", "w")
    r1_filenames_dict[i] = r1_filename
    r2_filenames_dict[i] = r2_filename

#Remember to convert to gz.open for actual read
with gz.open(r1, "rt") as R1, gz.open(r2, "rt") as R2, gz.open(r3, "rt") as R3, gz.open(r4, "rt") as R4, open("r1_unknown", "wt") as r1unknown, open("r2_unknown", "wt") as r2unknown, open("r1_hopped", "wt") as r1hopped, open("r2_hopped", "wt") as r2hopped:
#gzip.open(r1, mode = 'rt') as R1, gzip.open(r2, mode = 'rt') as R2, gzip.open(r3, mode = 'rt') as R3, gzip.open(r4, mode = 'rt'):
    #Create a loop to run through my files and collect one record at a time from each
    while True:
        #Read1 record
        header = R1.readline().strip()
        #In the case that we reach the end of the file
        if header == "":
            break
        sequence = R1.readline().strip()
        plus = R1.readline().strip()
        qscores = R1.readline().strip()
        r1_record = (header, sequence, plus, qscores)

        #Read2 record
        header = R2.readline().strip()
        sequence = R2.readline().strip()
        plus = R2.readline().strip()
        qscores = R2.readline().strip()
        r2_record = (header, sequence, plus, qscores)

        #Read3 record
        header = R3.readline().strip()
        sequence = R3.readline().strip()
        plus = R3.readline().strip()
        qscores = R3.readline().strip()
        r3_record = (header, sequence, plus, qscores)

        #Read4 record
        header = R4.readline().strip()
        sequence = R4.readline().strip()
        plus = R4.readline().strip()
        qscores = R4.readline().strip()
        r4_record = (header, sequence, plus, qscores)
        #print(r1_record, r2_record, r3_record, r4_record)

        #Check each record 
        #if the indexes are in the list of indexes
        if r2_record[1] in index_list and reverse_complement(r3_record[1]) in index_list: #problem must be either here or 
            #print(r2_record[1], "and",  reverse_complement(r3_record[1]), "FOUND")
            #print(index_quality_check(r2_record[1], r2_record[3]), index_quality_check(r3_record[1], r3_record[3]))
        #if the index meets the minimum quality check
            if index_quality_check(r2_record[1], r2_record[3]) == True and index_quality_check(r3_record[1], r3_record[3]) == True: #problem must be here
                #print("PASSED Q CHECK")
                #now check to see if the indexes match or are hopped
                #if they are matched write to the matched file
                if matched_hopped(r2_record[1], r3_record[1]) == True:
                    #print("MATCH", r2_record[1],reverse_complement(r3_record[1]))
                    #make sure to count this in the dictionary
                    index_pair_instances[r2_record[1], reverse_complement(r3_record[1])] += 1
                    r1_filenames_dict[i].write(r1_record[0] +" " + r2_record[1] + reverse_complement(r3_record[1]) + "\n" + r1_record[1] + "\n" + r1_record[2] + "\n" + r1_record[3] + "\n")
                    r2_filenames_dict[i].write(r4_record[0] +" " + reverse_complement(r3_record[1]) + r2_record[1] + "\n" + r4_record[1] + "\n" + r4_record[2] + "\n" + r4_record[3] + "\n")
                    matched_count +=1
                    index_dict[r2_record[1]] +=1
                #else the indexes are hopped write to the hopped file
                else:
                    #increment the dictionary
                    index_pair_instances[r2_record[1], reverse_complement(r3_record[1])] += 1
                    r1hopped.write(r1_record[0] +" " + r2_record[1] + reverse_complement(r3_record[1]) + "\n" + r1_record[1] + "\n" + r1_record[2] + "\n" + r1_record[3] + "\n")
                    r2hopped.write(r4_record[0] +" " + reverse_complement(r3_record[1]) + r2_record[1] + "\n" + r4_record[1] + "\n" + r4_record[2] + "\n" + r4_record[3] + "\n")
                    hopped_count +=1
            #the indexes don't meet the minimum quality score
            else: 
                r1unknown.write(r1_record[0] + " " + r2_record[1] + reverse_complement(r3_record[1]) + "\n" + r1_record[1] + "\n" + r1_record[2] + "\n" + r1_record[3] + "\n")
                r2unknown.write(r4_record[0] +" " + reverse_complement(r3_record[1]) + r2_record[1] + "\n" + r4_record[1] + "\n" + r4_record[2] + "\n" + r4_record[3] + "\n")
                unknown_count +=1
        #else the indexes aren't in the known index list, write to the unknown file
        else:
            r1unknown.write(r1_record[0] +" " + r2_record[1] + reverse_complement(r3_record[1]) + "\n" + r1_record[1] + "\n" + r1_record[2] + "\n" + r1_record[3] + "\n")
            r2unknown.write(r4_record[0] +" " + reverse_complement(r3_record[1]) + r2_record[1] + "\n" + r4_record[1] + "\n" + r4_record[2] + "\n" + r4_record[3] + "\n")
            unknown_count += 1

#Sort the index_pair_instances dictionary prior to writing it to output
index_pair_instances_ordered = sorted(index_pair_instances, key = index_pair_instances.get, reverse = True)

#Give the user some helpful output
print("Index Cuttoff Score (Per Base): ", minimum_index_score)
print("Biological Cuttoff Score (Read Average): ", minimum_biological_read_average)
print("Number of matched reads: ", matched_count)
print("Number of hopped reads: ", hopped_count)
print("Number of unknown reads: ", + unknown_count)
print("Percent Matched: ", ((matched_count / (matched_count + hopped_count + unknown_count))*100))
print("Percent Hopped: ",  ((hopped_count / (matched_count + hopped_count + unknown_count))*100))
for i in index_dict:
    print("Percent ", i, ":", (index_dict[i]/matched_count)*100)   
for I in index_pair_instances:
        print(I, ":", index_pair_instances[I])

#Store the above in an output file for reference later
with open("Demultiplex_Statistics.md","w") as output:
    output.write("Index Cuttoff Score (Per Base): " + str(minimum_index_score) + "\n")
    output.write("Biological Cuttoff Score (Read Average): " + str(minimum_biological_read_average) + "\n")
    output.write("Number of matched reads: " + str(matched_count) + "\n")
    output.write("Number of hopped reads: " + str(hopped_count) + "\n")
    output.write("Number of unknown reads: " + str(unknown_count) + "\n")
    output.write("Percent Matched: " + str(((matched_count / (matched_count + hopped_count + unknown_count))*100)) + "\n")
    output.write("Percent Hopped: " + str(((hopped_count / (matched_count + hopped_count + unknown_count))*100)) + "\n")
    for i in index_dict:
        output.write("Percent " + str(i) + " : " + str((index_dict[i]/matched_count)*100) + "\n")
    for I in index_pair_instances:
        output.write(str(I) + " : " + str(index_pair_instances[I]) + "\n")











# index2 = "ATCG"
# index3 = "CGAT"

# if matched_hopped(index2, index3) == True:
#     print("Matched")
# elif matched_hopped(index2, index3) == False:
#     print("Hopped")
# else:
#     print("We got issues")