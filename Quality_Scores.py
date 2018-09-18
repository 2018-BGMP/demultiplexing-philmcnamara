#!/usr/bin/env python3

import numpy as np
import matplotlib as mpl
#Necessary import to save plots on Talapas
mpl.use('Agg')
import matplotlib.pyplot as plt
import gzip
import argparse

#Argparse functionality
def get_arguments():
    parser = argparse.ArgumentParser(description = "Quality Score Distribution Visualizer")
    parser.add_argument("-l", "--read_length", help="Length of the Illumina Reads", required=True, type=int)
    parser.add_argument("-r", "--reads", help="Number of reads in input file", required=True, type = int)
    parser.add_argument("-f", "--file", help="Input file (FASTQ Format)", required=True, type = str)
    return parser.parse_args()
args = get_arguments()

read_length = args.read_length
file = args.file
lines = args.reads

#Phred-33 encoding conversion from ASCII to QScore
def convert_phred(letter):
    """Converts a single character into a phred score"""
    return ord(letter) - 33

#Initialize arrays, use running sum to save memory
sum_qscores = np.zeros(read_length, dtype=float)
mean_qscore = np.zeros(read_length, dtype=float)

#Line Counter
LC = 1

#Put the quality scores in the numpy array
with gzip.open(file, "rt") as f:
    for line in f:
        line = line.strip()
        if(LC%4 == 0):
            for x in range(read_length):
                sum_qscores[x] += convert_phred(line[x])
        if(LC % 100000 == 0):
            print("Working on line: " + str(LC))
        LC += 1

#Populate the array of mean qscores
for position in range(read_length):
    mean_qscore[position] = (sum_qscores[position] / lines)

#Rename file name, slicing is not ideal
file = file[49:55]

#Print the final array
print(file)
print(mean_qscore)

#Save and print figures
fig=plt.figure(figsize=(18, 16))
plt.rcParams.update({'font.size': 22})
plt.xlabel("Base Position")
plt.ylabel("Phred Quality Score")
plt.title("Mean Phred Quality Scores for " + file)
plt.plot(range(read_length), mean_qscore)
plt.savefig(file + "_Mean.png")
