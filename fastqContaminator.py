import argparse
import sys
import subprocess

# Purpose:
# Simulates a fastq file that has been produced from cross-contaminated NGS samples.  The user provides
# paired-end fastq sample files and paired-end fastq 'contaminating' files, specifies the proportion
# of contamination required - the script randomly selects reads form the 'contaminating' file and
# sample files in the correct proportions and concatenates them into a single file

# Requirements:
# Uses the tool fastq-sample to sample random reads from a fastq file
# (See https://homes.cs.washington.edu/~dcjones/fastq-tools/fastq-sample.html).
# fastq-tools-0.8 requires the libraries zlib (https://zlib.net/) and pcre (https://www.pcre.org/) to be installed.
# See the readme.txt for instructions on installing these C libraries and fastq-tools-0.8

# Usage:
# fastqContaminator.py -f test_WES21_Pan493_S9_R1_001.fastq test_WES21_Pan493_S9_R2_001.fastq -c test_WES21_Pan493_S10_R2_001.fastq test_WES21_Pan493_S10_R2_001.fastq -p 0.5 -s 5

###################### Argument Parser ######################

# Parses command line for arguments
# Usage:
# fastqContaminator.py -f test_WES21_Pan493_S9_R1_001.fastq test_WES21_Pan493_S9_R2_001.fastq -c test_WES21_Pan493_S10_R2_001.fastq test_WES21_Pan493_S10_R2_001.fastq -p 0.5 -s 5

parser = argparse.ArgumentParser(
    description='Uses the C program fastq-sample from the fastq-tools suite to simulate a fastq file with x% contamination from another sample.')
parser.add_argument('-sampleFiles', '-f', type=str, nargs=2,
                    help='Two paired-end fastq files containing the sample data')
parser.add_argument('-contaminantFiles', '-c', type=str, nargs=2,
                    help='Two paired-end fastq files containing the  contaminating data')
parser.add_argument('-proportionContaminant', '-p', type=float, nargs='*',
                    help='Array specifying the proportion (>0 & <=0.5) of contamination required.')
parser.add_argument('-seed', '-s', type=int,
                    help='Seed the random number generator. The same seed on the same data set will produce the same random sample.')
args = parser.parse_args()

###################### Functions ######################

# Function counts the lines in a fastq file:
def countFileLines(myFile):
    num_lines = sum(1 for line in open(myFile))
    return num_lines


# Check the input files are correct:
def checkInputFiles(sampleReads1, sampleReads2, contaminantReads1, contaminantReads2):
    # Check that both the supplied sample files contain the same number of paired-end reads
    if (countFileLines(sampleReads1) != countFileLines(sampleReads2)):
        print "The 2 paired-end sample files provided have a different number of reads/lines."
        sys.exit()
    # Check that both supplied contaminant files contain the same number of paired-end reads
    if (countFileLines(contaminantReads1) != countFileLines(contaminantReads2)):
        print "The 2 paired-end contamination files provided have different number of reads/lines."
        sys.exit()


# Function calculates the number of lines required from the sample/contaminating file to simulate x% contamination when added to sample file.
def calculateLinesRequired(sampleReads1, contaminantReads1, proportionContaminant):
    # Calculate the number of lines in the current sample file:
    numOfLines = countFileLines(sampleReads1)

    # Calculate the number of lines to sample from sample file to add to new file:
    numOfSampleLines = int(numOfLines * (1 - proportionContaminant))
    numOfSampleReads = int(numOfSampleLines / 4)

    # Calculate the number of lines to sample from contaminating file to equal required contamination:
    numOfContaminantLines = int(numOfLines * proportionContaminant)
    numOfContaminantReads = int(numOfContaminantLines / 4)

    # Check that contaminating file is sufficiently large:
    if (countFileLines(contaminantReads1) < numOfContaminantLines):
        print "Contamination file too small to provide the %d lines required for %f contamination." % (
        numOfContaminantLines, proportionContaminant)
        sys.exit()
    return numOfLines, numOfSampleReads, numOfContaminantReads


def sampleFastqReads(file1, file2, seedNum, numOfLinesRequired, myPrefix):
    myCommand = "fastq-sample -n %d -o %s -s %d %s %s" % (numOfLinesRequired, myPrefix, seedNum, file1, file2)
    return_code = subprocess.call(myCommand, shell=True)
    return return_code


def simulateContamination(sampleReads1, sampleReads2, contaminantReads1, contaminantReads2,
                          proportionContaminant, numOfSampleReads, numOfContaminantReads, seedNum):
    # Create prefix names based on % contamination:
    temp_sample = "temp_sample%s" % int((proportionContaminant * 100))
    temp_contaminant = "temp_contaminant%s" % int((proportionContaminant * 100))

    # Get reads from the sample files:
    sampleFastqReads(sampleReads1, sampleReads2, seedNum, numOfSampleReads, temp_sample)
    
    # Get reads from the contaminant files:
    sampleFastqReads(contaminantReads1, contaminantReads2, seedNum, numOfContaminantReads, temp_contaminant)

    # Combine the two produced files into one file and remove temporary files:
    myCommand = "cat %s.1.fastq %s.1.fastq > conPercent%s_%s &&" \
                " cat %s.2.fastq %s.2.fastq > conPercent%s_%s &&" \
                " rm %s.1.fastq %s.1.fastq %s.2.fastq %s.2.fastq" \
                % (temp_sample, temp_contaminant, int((proportionContaminant * 100)), sampleReads1,
                   temp_sample, temp_contaminant, int((proportionContaminant * 100)), sampleReads2,
                   temp_sample, temp_contaminant, temp_sample, temp_contaminant)
    return_code = subprocess.call(myCommand, shell=True)
    return return_code

###################### Main ####################

# Extract arguments from the command line:
# Filenames of each of the two sample fastq files containing the paired-end reads.
sampleReads1 = args.sampleFiles[0]
sampleReads2 = args.sampleFiles[1]

# Filenames of each of the two contaminating fastq files containing the paired-end reads.
contaminantReads1 = args.contaminantFiles[0]
contaminantReads2 = args.contaminantFiles[1]

proportionContaminant = args.proportionContaminant
seedNum = args.seed

# Sanity check on input
checkInputFiles(sampleReads1, sampleReads2, contaminantReads1, contaminantReads2)

#Reiterate over supplied array of user requsted proportions
for prop in proportionContaminant:

    # Calculate number of reads required from each file
    numOfLines, numOfSampleReads, numOfContaminantReads = calculateLinesRequired(sampleReads1, contaminantReads1,
                                                                             prop)

    # Simulate a file with x% contaminantion:
    simulateContamination(sampleReads1, sampleReads2, contaminantReads1, contaminantReads2,
                      prop, numOfSampleReads, numOfContaminantReads, seedNum)
