#!/usr/bin/env python3


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
import argparse

# Arguments
parser = argparse.ArgumentParser(description="Divide whole genomes into the genic regions")
parser.add_argument('-i', '--inputFile', help='File containing aligned genomes in fasta format', required=True)
parser.add_argument('-o', '--outputFilePrefix', help='Prefix for all gene.out files', required=True)
parser.add_argument('-g', '--geneGuidelines',
                    help='CSV file that contains geneic regions and intergenic regions labeled. EX: inter,21556,21561. inter for intergenic, start pos, end pos',
                    required=True)
args = parser.parse_args()

wholeGenFile = args.inputFile
outputPrefix = args.outputFilePrefix


def findStops(seq, regionStart, regionEnd):
    stop1 = "TAA"
    stop2 = "TAG"
    stop3 = "TGA"
    for i in range(regionStart, regionEnd, 3):
        codon = seq[i:i + 3]
        if (stop1 in codon) or (stop2 in codon) or (stop3 in codon):
            return False  # Stop codon found
    return True  # No stop codon found


geneGuidelines = []
with open(args.geneGuidelines, newline='') as file:
    reader = csv.reader(file)
    for row in reader:
        geneGuidelines.append(row)

badSeqs = []
for sublist in geneGuidelines:
    if sublist[0] != "inter":
        for record in SeqIO.parse(wholeGenFile, "fasta"):
            seq = str(record.seq)
            regionStart = int(sublist[1])
            regionEnd = int(sublist[2]) - 3
            if not findStops(seq, regionStart, regionEnd):
                badSeqs.append(record.id)
                print(record.id, "is a bad seq")


for sublist in geneGuidelines:
    length = 60
    if sublist[0] != "inter":
        handle = outputPrefix + sublist[0]
        with open(handle, 'w') as outFile:
            for record in SeqIO.parse(wholeGenFile, "fasta"):
                seq = str(record.seq)
                regionStart = int(sublist[1])
                regionEnd = int(sublist[2])

                # print("Region: ", sublist[0])
                if record.id in badSeqs:
                    # print("Found", record.id, "in badseqs")
                    continue

                goi = seq[regionStart:regionEnd - 3]
                # goiRec = SeqRecord(Seq(goi), id=record.id, name=record.name)

                # print(f"Writing sequence: {record.id} to {sublist[0]} file")
                # write only the id and seq
                outFile.write(f">{record.id}\n")
                for i in range(0, len(goi), length):
                    outFile.write(goi[i:i + length] + '\n')

                # SeqIO.write(goiRec,outFile, "fasta")





