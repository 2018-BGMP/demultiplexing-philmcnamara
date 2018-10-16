#! /usr/bin/env python3

import os
import argparse
import statistics
import gzip


def get_arguments():
    parser = argparse.ArgumentParser(description="Demultiplexer")
    parser.add_argument("-i", "--index_cutoff",
                        help="minimum average quality score for an index \
                        (default 30)",
                        required=False, type=int, nargs="?", default=30)
    parser.add_argument("-r", "--read_cutoff",
                        help="minimum average quality score for a read \
                        (default 25)",
                        required=False, type=int, nargs="?", default=25)
    parser.add_argument("-f", "--files", help="Input the multiplexed FASTQ \
                        files in this order: Forward Read  Forward Index \
                        Reverse Read  Reverse Index",
                        required=True, type=str, nargs=4)
    return parser.parse_args()


args = get_arguments()

files = args.files
index_qscore_cutoff = args.index_cutoff
read_qscore_cutoff = args.read_cutoff

os.makedirs("demultiplexed")

# GLOBAL VARIABLES/DATA STRUCTURES

# I know this is not ideal but I couldn't figure out how
# to open and name files in a loop
B1_R1 = open("demultiplexed/B1_R1.fq", "w")
B1_R2 = open("demultiplexed/B1_R2.fq", "w")
A5_R1 = open("demultiplexed/A5_R1.fq", "w")
A5_R2 = open("demultiplexed/A5_R2.fq", "w")
C1_R1 = open("demultiplexed/C1_R1.fq", "w")
C1_R2 = open("demultiplexed/C1_R2.fq", "w")
B9_R1 = open("demultiplexed/B9_R1.fq", "w")
B9_R2 = open("demultiplexed/B9_R2.fq", "w")
C9_R1 = open("demultiplexed/C9_R1.fq", "w")
C9_R2 = open("demultiplexed/C9_R2.fq", "w")
C3_R1 = open("demultiplexed/C3_R1.fq", "w")
C3_R2 = open("demultiplexed/C3_R2.fq", "w")
B3_R1 = open("demultiplexed/B3_R1.fq", "w")
B3_R2 = open("demultiplexed/B3_R2.fq", "w")
C4_R1 = open("demultiplexed/C4_R1.fq", "w")
C4_R2 = open("demultiplexed/C4_R2.fq", "w")
A11_R1 = open("demultiplexed/A11_R1.fq", "w")
A11_R2 = open("demultiplexed/A11_R2.fq", "w")
C7_R1 = open("demultiplexed/C7_R1.fq", "w")
C7_R2 = open("demultiplexed/C7_R2.fq", "w")
B2_R1 = open("demultiplexed/B2_R1.fq", "w")
B2_R2 = open("demultiplexed/B2_R2.fq", "w")
A1_R1 = open("demultiplexed/A1_R1.fq", "w")
A1_R2 = open("demultiplexed/A1_R2.fq", "w")
B7_R1 = open("demultiplexed/B7_R1.fq", "w")
B7_R2 = open("demultiplexed/B7_R2.fq", "w")
A3_R1 = open("demultiplexed/A3_R1.fq", "w")
A3_R2 = open("demultiplexed/A3_R2.fq", "w")
B4_R1 = open("demultiplexed/B4_R1.fq", "w")
B4_R2 = open("demultiplexed/B4_R2.fq", "w")
A12_R1 = open("demultiplexed/A12_R1.fq", "w")
A12_R2 = open("demultiplexed/A12_R2.fq", "w")
C10_R1 = open("demultiplexed/C10_R1.fq", "w")
C10_R2 = open("demultiplexed/C10_R2.fq", "w")
A2_R1 = open("demultiplexed/A2_R1.fq", "w")
A2_R2 = open("demultiplexed/A2_R2.fq", "w")
C2_R1 = open("demultiplexed/C2_R1.fq", "w")
C2_R2 = open("demultiplexed/C2_R2.fq", "w")
A10_R1 = open("demultiplexed/A10_R1.fq", "w")
A10_R2 = open("demultiplexed/A10_R2.fq", "w")
B8_R1 = open("demultiplexed/B8_R1.fq", "w")
B8_R2 = open("demultiplexed/B8_R2.fq", "w")
A7_R1 = open("demultiplexed/A7_R1.fq", "w")
A7_R2 = open("demultiplexed/A7_R2.fq", "w")
B10_R1 = open("demultiplexed/B10_R1.fq", "w")
B10_R2 = open("demultiplexed/B10_R2.fq", "w")
A8_R1 = open("demultiplexed/A8_R1.fq", "w")
A8_R2 = open("demultiplexed/A8_R2.fq", "w")
unk_r1 = open("demultiplexed/unknown_read1.fq", "w")
unk_r2 = open("demultiplexed/unknown_read2.fq", "w")

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

out_files = {
    "B1": [B1_R1, B1_R2],
    "A5": [A5_R1, A5_R2],
    "C1": [C1_R1, C1_R2],
    "B9": [B9_R1, B9_R2],
    "C9": [C9_R1, C9_R2],
    "C3": [C3_R1, C3_R2],
    "B3": [B3_R1, B3_R2],
    "C4": [C4_R1, C4_R2],
    "A11": [A11_R1, A11_R2],
    "C7": [C7_R1, C7_R2],
    "B2": [B2_R1, B2_R2],
    "A1": [A1_R1, A1_R2],
    "B7": [B7_R1, B7_R2],
    "A3": [A3_R1, A3_R2],
    "B4": [B4_R1, B4_R2],
    "A12": [A12_R1, A12_R2],
    "C10": [C10_R1, C10_R2],
    "A2": [A2_R1, A2_R2],
    "C2": [C2_R1, C2_R2],
    "A10": [A10_R1, A10_R2],
    "B8": [B8_R1, B8_R2],
    "A7": [A7_R1, A7_R2],
    "B10": [B10_R1, B10_R2],
    "A8": [A8_R1, A8_R2]
}

# Base-pairing for finding reverse complement
base_pairing = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

# Properly Matched Tracker Dictionary
properly_matched = dict()
for value in index_sequences.values():
    properly_matched[value] = 0

# Number of reads that pass quality score cutoff
quality_reads = 0
# Number of reads that pass quality score AND have correct indexes
# Reads that pass below AND aren't index-hopped will be written to output and
# counted in the "properly-matched" dictionary
good_indexes = 0

# Total number of sequences processed
total_sequence_count = 0


# FUNCTIONS


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

# FILE PROCESSING BODY


print("Starting demultiplexing with minimum average read quality score of "
      + str(read_qscore_cutoff) +
      " and minimum average index quality score of "
      + str(index_qscore_cutoff))

with gzip.open(files[0], "rt") as Read1_file, \
        gzip.open(files[1], "rt") as Index1_file, \
        gzip.open(files[2], "rt") as Read2_file, \
        gzip.open(files[3], "rt") as Index2_file:
    for line in Read1_file:
        if(total_sequence_count % 10000000 == 0):
            print("Now processing sequence: " + str(total_sequence_count))
        Read1 = []
        Index1 = []
        Read2 = []
        Index2 = []
        Read1.append(line.strip())
        for i in range(3):
            line = Read1_file.readline()
            Read1.append(line.strip())
        for i in range(4):
            line = Index1_file.readline()
            Index1.append(line.strip())
        for i in range(4):
            line = Read2_file.readline()
            Read2.append(line.strip())
        for i in range(4):
            line = Index2_file.readline()
            Index2.append(line.strip())
        total_sequence_count += 2
        # Discard reads if they are below a quality score cutoff
        if((mean_read_quality(Read1[3])) > read_qscore_cutoff
           and (mean_read_quality(Read2[3]) > read_qscore_cutoff
                and (mean_read_quality(Index1[3]) > index_qscore_cutoff)
                and (mean_read_quality(Index2[3]) > index_qscore_cutoff))):
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
                    for i in range(4):
                        out_files[index_name][0].write(Read1[i])
                        # Write the index sequence in the header
                        if(i == 0):
                            out_files[index_name][0].write(":" + Index1[1])
                        out_files[index_name][0].write('\n')
                    for i in range(4):
                        out_files[index_name][1].write(Read2[i])
                        if(i == 0):
                            out_files[index_name][0].write(":" + Index1[1])
                        out_files[index_name][1].write('\n')
                else:
                    for i in range(4):
                        unk_r1.write(Read1[i])
                        unk_r1.write('\n')
                        unk_r2.write(Read2[i])
                        unk_r2.write('\n')
            else:
                for i in range(4):
                    unk_r1.write(Read1[i])
                    unk_r1.write('\n')
                    unk_r2.write(Read2[i])
                    unk_r2.write('\n')

unk_r1.close()
unk_r2.close()

# RESULTS REPORTING
matched_reads = 0
for key, value in properly_matched.items():
    matched_reads += value
with open("results_testing.txt", "w") as o:
    o.write("Total Reads Processed: " + str(total_sequence_count) + "\n\n")
    low_qual = total_sequence_count - quality_reads
    percent_low_qual = round((low_qual/total_sequence_count) * 100, 2)
    o.write("Reads discarded due to low quality: "
            + str(low_qual) + " (" + str(percent_low_qual) + "%)\n\n")
    unknown = quality_reads - good_indexes
    percent_unknown = round((unknown/quality_reads) * 100, 2)
    o.write("Reads with unknown indexes: " + str(unknown) + " (" +
            str(percent_unknown) + "%)\n\n")
    o.write("Demultiplexed read counts by index\n")
    for key, value in properly_matched.items():
        percent_by_index = round((value/matched_reads) * 100, 2)
        o.write(key + ": " + str(value) + " sequences  "
                + str(percent_by_index) + "%\n")
    o.write("\n")
    o.write("Index hopping is only assessed for reads above the quality score "
            "filter and with two correct indexes\n")
    hopping = round((1 - (matched_reads/good_indexes)) * 100, 4)
    o.write("Percentage Index Hopping: " + str(hopping) + "%\n")
