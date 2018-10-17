# Demultiplexing

All initial work from Part 1 is in its own repository

## Unit Tests

*unit_tests* folder holds test files for 8 cases, running the script on these cases gives **results_unit_test.txt**

Output from running demultiplexing.py on unit tests is stored in *demultiplexed* directory

All .fq files created are blank except A1 (the index in my test cases) and the two unknown files

## 2017 Sequencing Data

Running the script on the 2017 sequencing data gives **results.txt**, demultiplexed .fq files are stored on Talapas

Runtime was 28 hours 32 minutes on 8 cores

## Update - Code review

After code review/cleanup my runtime was 2 hours 30 minutes on 1 core. The factor slowing me down the most was calculating the mean quality score of every read in addition the indexes. It was also slow because I was storing the data in an array and taking the mean. 

I made an update to search for the index in the dictionary instead of dictionary.keys(), which apparently does a linear search through a list instead of a hashed lookup.