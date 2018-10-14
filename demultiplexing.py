#!/usr/bin/env python3

import statistics
import os

# TODO Add argparse functionality (Include minimum quality score cutoff as optional parameters)

# Dictionary of Index Sequences
index_sequences = {
    "GTAGCGTA": "B1",
    "CGATCGAT": "A5",
    "GATCAAGG": "C1",
    "AACAGCGA": "B9",
    "TAGCCATG": "C9",
    "CGGTAATC": "C3",
    "CTCTGGAT": "B3",
    "TACCGGAT": "C4",
    "CTAGCTCA": "A11",
    "CACTTCAC": "C7",
    "GCTACTCT": "B2",
    "ACGATCAG": "A1",
    "TATGGCAC": "B7",
    "TGTTCCGT": "A3",
    "GTCCTAAG": "B4",
    "TCGACAAG": "A12",
    "TCTTCGAC": "C10",
    "ATCATGCG": "A2",
    "ATCGTGGT": "C2",
    "TCGAGAGT": "A10",
    "TCGGATTC": "B8",
    "GATCTTGC": "A7",
    "AGAGTCCA": "B10",
    "AGGATAGC": "A8"
}

# Properly Matched Tracker Dictionary
properly_matched = dict()
for value in index_sequences.values():
    properly_matched[value] = 0

# Number of reads that pass quality score cutoff
quality_reads = 0
# Number of reads that pass quality score AND have correct indexes
good_indexes = 0
# Reads that pass above AND aren't index-hopped will be written to output and
# counted in the "properly-matched" dictionary


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


# os.makedirs("output")

with open("test_files/miniR1.fq", "r") as R1file, \
        open("test_files/miniR2.fq", "r") as R2file, \
        open("test_files/miniR3.fq", "r") as R3file, \
        open("test_files/miniR4.fq", "r") as R4file:
    # for index in index_sequences.values():
    #     open("output/" + index + "_read1.fq", "w")
    #     open("output/" + index + "_read2.fq", "w")
    # open("output/unknown_read1.fq", "w")
    # open("output/unknown_read2.fq", "w")
    for line in R1file:
        Read1 = []
        Index1 = []
        Index2 = []
        Read2 = []
        Read1.append(line.strip())
        for i in range(3):
            line = R1file.readline()
            Read1.append(line.strip())
        for i in range(4):
            line = R2file.readline()
            Index1.append(line.strip())
        for i in range(4):
            line = R3file.readline()
            Index2.append(line.strip())
        for i in range(4):
            line = R4file.readline()
            Read2.append(line.strip())

        # Discard reads if they are below a quality score cutoff
        if((mean_read_quality(Read1[3])) > 25
           and (mean_read_quality(Read2[3]) > 25)
           and (mean_read_quality(Index1[3]) > 30)
           and (mean_read_quality(Index2[3]) > 30)):
            quality_reads += 2
            # Check if the indexes of the reads are accurate
            if(Index1[1] in index_sequences.keys()
                    and Index2[1] in index_sequences.keys()):
                good_indexes += 2
                # Check for index-hopping
                if(Index1[1] == Index2[1]):
                    index_name = index_sequences[Index1[1]]
                    properly_matched[index_name] += 2
                    print("Match")
                    # Write to the appropriate index-read file
                else:
                    pass
                    # Write to unknown index sequence files
                    print("Index-Hop")
            else:
                pass
                # Write to unknown index sequence files
                print("Bad Index")
        else:
            pass
            print("Low-Quality")
            # Results output:
            # TODO Format Results
