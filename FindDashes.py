from Bio import SeqIO
import argparse
import csv

parser = argparse.ArgumentParser(description= "Ensure aligned genomes are in the same ORF ")
parser.add_argument('-i', '--inputFile', help='Aligned genomes in fasta format')
parser.add_argument('-o', '--outputFile', help = 'File where valid genomes will be copied')
args = parser.parse_args()
fastaFile =  args.inputFile
printFile = args.outputFile


def copySeqs(inputFile, outputFile):
    numValid = 0
    numInvalid = 0
    with open(inputFile, 'r') as inputHandle, open(outputFile, 'w') as outputHandle:
        for record in SeqIO.parse(inputHandle, "fasta"):
            if isValid(record.seq):
                numValid += 1
                #print("hello", file=outputHandle)
                SeqIO.write(record, outputHandle, "fasta")
            else:
                numInvalid += 1
    print(f"There were {numValid} sequences printed. \nThere were {numInvalid} sequences rejected")

def isValid(seq):
    consecutive = 0
    keepSeq = True
    for char in seq:
        if char == '-':
            consecutive += 1
        else:
            if consecutive % 3 != 0:
                consecutive = 0
                keepSeq = False
            else:
                consecutive = 0
    return keepSeq


copySeqs(fastaFile, printFile)


