#! /usr/bin/env python3

import os
import argparse
import gzip

# ARGPARSE


def get_arguments():
    parser = argparse.ArgumentParser(description="Demultiplexer")
    parser.add_argument("-i", "--index_cutoff",
                        help="minimum average quality score for an index \
                        (default 30)",
                        required=False, type=int, nargs="?", default=30)
    parser.add_argument("-f", "--files", help="Input the multiplexed FASTQ \
                        files in this order: Forward Read  Forward Index \
                        Reverse Read  Reverse Index",
                        required=True, type=str, nargs=4)
    return parser.parse_args()


args = get_arguments()

files = args.files
index_qscore_cutoff = args.index_cutoff

os.makedirs("demultiplexed")

# GLOBAL VARIABLES/DATA STRUCTURES

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

# Dictionary to hold output file names as keys, file handles as values
out_files = {}

# For each index name, open a R1 and R2 file.
# This avoids opening/closing the files for each read/write operation
for value in index_sequences.values():
    out_files[value + "_R1"] = open("demultiplexed/" + value + "_R1.fq", "w")
    out_files[value + "_R2"] = open("demultiplexed/" + value + "_R2.fq", "w")

unk_r1 = open("demultiplexed/unknown_read1.fq", "w")
unk_r2 = open("demultiplexed/unknown_read2.fq", "w")

# Base-pairing for finding reverse complement, N is unchanged
base_pairing = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

# Dictionary to hold number of reads for each index
# that are high quality and not index hopped
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


def mean_index_quality(qscore_line):
    '''Coverts a string of quality score ASCII characters into a list,
    then calls convert_phred to create a second list of quality score integers.
    Calculates the mean value of the list and returns it.'''
    sum = 0
    for score in qscore_line:
        # Convert from Phred-33
        sum += ord(score) - 33
    return sum/8


def reverse_complement(sequence):
    '''Takes a DNA sequence string and returns the
    reverse complementary sequence'''
    reversed_sequence = sequence[::-1]
    reverse_complement = ""
    for base in reversed_sequence:
        reverse_complement += base_pairing[base]
    return reverse_complement

# FILE PROCESSING BODY


print("Starting demultiplexing with minimum average index quality score of "
      + str(index_qscore_cutoff))

with gzip.open(files[0], "rt") as Read1_file, \
        gzip.open(files[1], "rt") as Index1_file, \
        gzip.open(files[2], "rt") as Read2_file, \
        gzip.open(files[3], "rt") as Index2_file:
    # Progress indicator
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
        # Increment by 2 for forward and reverse, same for other counters below
        total_sequence_count += 2
        # Discard reads if they are below a quality score cutoff
        if((mean_index_quality(Index1[3]) > index_qscore_cutoff)
                and (mean_index_quality(Index2[3]) > index_qscore_cutoff)):
            quality_reads += 2
            # Use the reverse complement of Index 2 for comparison
            Index2_RC = reverse_complement(Index2[1])
            # Check if the indexes of the reads are accurate
            if(Index1[1] in index_sequences and Index2_RC in index_sequences):
                good_indexes += 2
                # Check for index-hopping
                if(Index1[1] == Index2_RC):
                    index_name = index_sequences[Index1[1]]
                    properly_matched[index_name] += 2
                    # Write each line of the read to output file
                    for i in range(4):
                        out_files[index_name + "_R1"].write(Read1[i])
                        # Write the index sequence in the @ header
                        if(i == 0):
                            out_files[index_name + "_R1"].write(":" + Index1[1])
                        out_files[index_name + "_R1"].write('\n')
                    for i in range(4):
                        out_files[index_name + "_R2"].write(Read2[i])
                        if(i == 0):
                            out_files[index_name + "_R2"].write(":" + Index1[1])
                        out_files[index_name + "_R2"].write('\n')
                # Index-hopping alternative
                else:

                    for i in range(4):
                        unk_r1.write(Read1[i] + "\n")
                        unk_r2.write(Read2[i] + "\n")
            # Bad index alternative
            else:
                for i in range(4):
                    unk_r1.write(Read1[i] + "\n")
                    unk_r2.write(Read2[i] + "\n")
        # Low index quality alternative
        else:
            for i in range(4):
                unk_r1.write(Read1[i] + "\n")
                unk_r2.write(Read2[i] + "\n")

# Close all files not automatically closed by "with"

unk_r1.close()
unk_r2.close()

for value in out_files.values():
    value.close()

# RESULTS REPORTING

# Sum the total number of reads with no index hopping
matched_reads = 0
for key, value in properly_matched.items():
    matched_reads += value
with open("results.txt", "w") as o:
    o.write("Total Reads Processed: " + str(total_sequence_count) + "\n\n")
    low_qual = total_sequence_count - quality_reads
    percent_low_qual = round((low_qual/total_sequence_count) * 100, 2)
    o.write("Reads with low quality score indexes: "
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
