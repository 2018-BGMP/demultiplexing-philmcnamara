#!/usr/bin/env python3

import statistics
import os
import argparse

# TODO Add argparse functionality (Include minimum quality score cutoff as optional parameters)


# def get_arguments():
#     parser = argparse.ArgumentParser(description="Demultiplexer")
#     parser.add_argument("-c", "--coverage_limit",
#                         help="k-mer coverage limit", required=True, type=int)
#     parser.add_argument("-f", "--file", help="Input file (FASTQ Format)", required=True, type=str)
#     return parser.parse_args()
#
#
# args = get_arguments()

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

# Base-pairing for finding reverse complement
base_pairing = {"A": "T", "T": "A", "G": "C", "C": "G"}

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

total_sequence_count = 0


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


def reverse_complement(sequence):
    '''Takes a DNA sequence string and returns the
    reverse complementary sequence'''
    reversed_sequence = sequence[::-1]
    reverse_complement = ""
    for base in reversed_sequence:
        reverse_complement += base_pairing[base]
    return reverse_complement


os.makedirs("output")

with open("test_files/miniR1.fq", "r") as R1file, \
        open("test_files/miniR2.fq", "r") as R2file, \
        open("test_files/miniR3.fq", "r") as R3file, \
        open("test_files/miniR4.fq", "r") as R4file:
    unk_r1 = open("output/unknown_read1.fq", "w")
    unk_r2 = open("output/unknown_read2.fq", "w")
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
        total_sequence_count += 2
        # Discard reads if they are below a quality score cutoff
        if((mean_read_quality(Read1[3])) > 25
           and (mean_read_quality(Read2[3]) > 25)
           and (mean_read_quality(Index1[3]) > 30)
           and (mean_read_quality(Index2[3]) > 30)):
            quality_reads += 2
            # Use the reverse complement of Index 2 for comparison
            Index2_RC = reverse_complement(Index2[1])
            # Check if the indexes of the reads are accurate
            if(Index1[1] in index_sequences.keys()
                    and Index2_RC in index_sequences.keys()):
                good_indexes += 2
                # Check for index-hopping
                if(Index1[1] == Index2_RC):
                    index_name = index_sequences[Index1[1]]
                    properly_matched[index_name] += 2
                    print("Match")
                    with open("output/{}_R1".format(index_name), "w") as f:
                        for i in range(4):
                            f.write(Read1[i])
                            f.write('\n')
                    with open("output/{}_R2".format(index_name), "w") as f:
                        for i in range(4):
                            f.write(Read2[i])
                            f.write('\n')
                else:
                    for i in range(4):
                        unk_r1.write(Read1[i])
                        unk_r1.write('\n')
                        unk_r2.write(Read2[i])
                        unk_r2.write('\n')
                    print("Index-Hop")
            else:
                for i in range(4):
                    unk_r1.write(Read1[i])
                    unk_r1.write('\n')
                    unk_r2.write(Read2[i])
                    unk_r2.write('\n')

unk_r1.close()
unk_r2.close()

# Results output:
# TODO Format Results

matched_reads = 0
for key, value in properly_matched.items():
    matched_reads += value
with open("results.txt", "w") as o:
    o.write("Total Reads Processed: " + str(total_sequence_count))
    o.write("Reads discarded due to low quality: "
            + str(total_sequence_count - quality_reads))
    o.write("Reads with unknown indexes: " + str(quality_reads - good_indexes))
    for key, value in properly_matched.items():
        o.write(key + ": " + str(value / matched_reads) + "%")
    o.write("Index hopping is only assessed for reads above the quality score \
            filter and with two correct indexes")
    o.write("Percentage Index Hopping: " +
            str(matched_reads / good_indexes) + "%")
