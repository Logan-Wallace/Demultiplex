#! usr/bin/env python

'''The purpose of this script is to parse a fastq file and return
a histogram of the quality score distribution'''

#Psuedocode
#import the bioinfo.py module
#call the init list function to create an empty list
#call the populate list function to loop through the fastqfile and 
#populate our list. 
#plot the result using matplotlib

#import necessary modules
import matplotlib.pyplot as plt
import bioinfo
import argparse

parser = argparse.ArgumentParser(description = "Arguments for demultiplex.py")
parser.add_argument("-f", "--fastq_filename", help = "The name of the fastqfile you wish to parse", type = str, default = "")
parser.add_argument("-r", "--read_length", help = "What is the length of the read within the fastq file we are reading from", type = int, default = 101)
parser.add_argument("-o", "--output_filename", help = "The name of the fastqfile you wish to parse", type = str, default = "output.png")
args = parser.parse_args()

#introduce my function from bioinfo.py
def init_list(lst: list, r: int, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    for n in range(r):
        lst.append(value)
    return lst

#Declare my variables
fastq = args.fastq_filename
r = args.read_length
output = args.output_filename
qscores: list = []
num_lines: int = 0

#Call my init list function
bioinfo.init_list(qscores)

#Call the populate_list function to fill the empty list with my qscores
qscores, num_lines = bioinfo.populate_list(fastq)

#print out the mean of qscores list just to make sure we are getting some actual readout
#print("# Base Pair" + "\t" + "Mean Quality Score" + "\n")
for i in range(len(qscores)):
        qscores[i] = (qscores[i] / (num_lines / 4)) 
#        print(str(i) + "\t" + str(qscores[i]))   

#Plot the results from our list as a histogram
plt.bar(range(len(qscores)), qscores)
plt.xlabel("Read Position")
plt.ylabel("Quality Score")
plt.title("Quality Score Distribution by Position in read")
plt.savefig(output + ".png")