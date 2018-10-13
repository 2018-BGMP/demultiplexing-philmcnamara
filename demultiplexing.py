#!/usr/bin/env python3

import statistics

# Dictionary of Index Sequences
index_sequences = {
    "B1": "GTAGCGTA",
    "A5": "CGATCGAT",
    "C1": "GATCAAGG",
    "B9": "AACAGCGA",
    "C9": "TAGCCATG",
    "C3": "CGGTAATC",
    "B3": "CTCTGGAT",
    "C4": "TACCGGAT",
    "A11": "CTAGCTCA",
    "C7": "CACTTCAC",
    "B2": "GCTACTCT",
    "A1": "ACGATCAG",
    "B7": "TATGGCAC",
    "A3": "TGTTCCGT",
    "B4": "GTCCTAAG",
    "A12": "TCGACAAG",
    "C10": "TCTTCGAC",
    "A2": "ATCATGCG",
    "C2": "ATCGTGGT",
    "A10": "TCGAGAGT",
    "B8": "TCGGATTC",
    "A7": "GATCTTGC",
    "B10": "AGAGTCCA",
    "A8": "AGGATAGC",
}

# Properly Matched Tracker Dictionary
properly_matched = dict()
for key in index_sequences.keys():
    properly_matched[key] = 0

# Number of reads that pass quality score cutoff
quality_reads = 0


def read_lines(file):
    '''Takes a FASTQ input file and reads four lines, stripping newline characters
    and saving them to strings prefixed by their file name:
    file_header, file_sequence, file_plus, file_quality'''
    lines = []
    # Read in the first 4 lines, 1 sequence
    for i in range(4):
        lines.append(file.readline())
    return lines


def convert_phred(char):
    '''Converts an ASCII phred score to a numeric quality score value'''
    return ord(char) - 33


def mean_read_quality(qscore_line):
    '''Coverts a string of quality score ASCII characters into a list,
    then calls convert_phred to create a second list of quality score integers.
    Calculates the mean value of the list and returns it.'''
    converted_list = []
    for score in qscore_line:
        score = convert_phred(score)
        converted_list.append(score)
    return statistics.mean(converted_list)


with open("test_files/R1.fq", "r") as R1, \
        open("test_files/R2.fq", "r") as R2, \
        open("test_files/R3.fq", "r") as R3, \
        open("test_files/R4.fq", "r") as R4:
    print(R1.readline())
    print(R2.readline())
    print(R3.readline())
    print(R4.readline())
