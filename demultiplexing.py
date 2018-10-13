#!/usr/bin/env python3

import statistics

# TODO Add argparse functionality
# Include minimum quality score cutoff as optional parameters

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


def read_sequence(file):
    '''Takes a FASTQ input file and reads four lines, stripping newline characters
    and saving them to strings prefixed by their file name:
    file_header, file_sequence, file_plus, file_quality'''
    lines = []
    # Read in the first 4 lines, 1 sequence
    for i in range(4):
        lines.append(file.readline().strip())
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
    # TODO Fix this to make it read the whole file without going over
    # TODO Open 50 files for writing
    Read1 = read_sequence(R1)
    Index1 = read_sequence(R2)
    Read2 = read_sequence(R4)
    Index2 = read_sequence(R3)
    print(Index1[1])
    print(Index2[1])

    # Discard reads if they are below a quality score cutoff
    if((mean_read_quality(Read1[3])) > 25 and (mean_read_quality(Read2[3]) > 25) and (mean_read_quality(Index1[3]) > 30) and (mean_read_quality(Index2[3]) > 30)):
        quality_reads += 2
        # Check if the indexes of the reads are accurate
        if(Index1[1] in index_sequences.keys() and Index2[1] in index_sequences.keys()):
            good_indexes += 2
            # Check for index-hopping
            if(Index1[1] == Index2[1]):
                index_name = index_sequences[Index1[1]]
                properly_matched[index_name] += 2
                # Write to the appropriate index-read file
            else:
                # Write to unknown index sequence files
                pass
        else:
            # Write to unknown index sequence files
            pass

# Results output:
