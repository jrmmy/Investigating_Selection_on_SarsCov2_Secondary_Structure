#!/usr/bin/env python3

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Align
from Bio import Entrez
import csv
import argparse

# Arguments
parser = argparse.ArgumentParser(description= "Ensure aligned genomes are in the same ORF, removes indels from non coding regions")
parser.add_argument('-i', '--inputFile', help='Aligned genomes in fasta format', required = True)
parser.add_argument('-o', '--outputFile', help = 'File where valid genomes will be copied', required = True)
parser.add_argument('-g', '--geneGuidelines', help = 'CSV file that contains geneic regions and intergenic regions labeled. EX: inter,21556,21561. inter for intergenic, start pos, end pos', required = True)
args = parser.parse_args()
fastaFile =  args.inputFile
printFile = args.outputFile

# Make list of lists for gene guidelines
geneGuidelines = []
with open(args.geneGuidelines, newline= '') as file:
    reader = csv.reader(file)
    for row in reader:
        geneGuidelines.append(row)

Entrez.email = 'jd2669@msstate.edu'
db = "nucleotide"
accession_number = "MN908947"

# Use Entrez.efetch to get handle from NCBI
handle = Entrez.efetch(db=db, id=accession_number, rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()
refGenome = record.seq

# Function to replace other two
def editAlignment(ref, wild, genes):
    offset = 0
    print("Editing Alignment", '\n')
    insertsLoc = []

    # Looking for insertions in genome of interest
    for sublist in genes:
        print(f"Working on gene: {sublist[0]} \n")
        regionStart = int(sublist[1])
        regionEnd = int(sublist[2])
        pos = ref.find("-", regionStart + offset, regionEnd + offset)
        while pos > 0:
            insertLength = 0
            #print('.find method worked correctly')
            while ref[pos] == "-":
                insertLength += 1
                insertsLoc.append(pos)
                print('Appending to insertsLoc \n')
                pos += 1
                offset += 1
            if sublist[0] != 'inter' and insertLength % 3 != 0:
                print("invalid sequence \n")
                print("Insertion location: ", insertsLoc, " - size = ", insertLength)
                return False, ""
            pos = ref.find("-", pos, regionEnd + offset)
    insertsLoc.sort(reverse=True)
    charWild = list(wild)
    #     #Loop through the locations marked 
    for i in insertsLoc:
        charWild.pop(i)
    wild = ''.join(charWild)
    print("All insertion locations:", insertsLoc)

    # Now looking for deletions in genome of interest
    for sublist in genes:
        if sublist[0] != "inter":
            regionStart = int(sublist[1])
            regionEnd = int(sublist[2])
            pos = wild.find("-", regionStart, regionEnd)
            while pos > 0:
                deletionLength = 0
                while wild[pos] == "-":
                    deletionLength += 1
                    pos += 1
                if deletionLength % 3 != 0:
                    print("Invalid Sequence \n")
                    print(f"Deletion size: {deletionLength}")
                    return False, ""
                pos = wild.find("-", pos, regionEnd)

    return True, str(wild)
                


aligner = Align.PairwiseAligner(mode='global', match_score=2, mismatch_score=-1, 
                          open_gap_score=-10, extend_gap_score=-0.5, 
                          target_end_gap_score = 0.0, query_end_gap_score = 0.0)

numValid = 0
numInvalid = 0

# Main loop
with open(fastaFile, 'r') as inputHandle, open(printFile, 'w') as outputHandle:
    print("Opened files successfully")
    # Loop through each genome in the input file
    for record in SeqIO.parse(inputHandle, "fasta"):

        print(f"Aligning sequence: {record}")
        print(' ')
        
        # Align w/ pairwiseAligner
        alignments = aligner.align(refGenome, record)

        print("Number of alignments returned: ", len(alignments))

        success, updatedSeq = editAlignment(alignments[0][0], alignments[0][1], geneGuidelines)
        # If find dashes returns true, write the sequence to the output file
        if success == True :
            print("Success is True")
            updatedRecord = SeqRecord(Seq(updatedSeq), id = record.id, name = record.name)    
            numValid += 1
            #print("hello", file=outputHandle)
            SeqIO.write(updatedRecord, outputHandle, "fasta")
        else:
            print("Success is False")
            numInvalid += 1
    print(f"There were {numValid} sequences printed. \nThere were {numInvalid} sequences rejected")