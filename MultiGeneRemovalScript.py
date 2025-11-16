#!/mnt/home/software/conda-envs/biopython/bin/python



from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import csv

def count_undef(seq):
    i = 0
    x = 0
    for i in range(len(seq)):
        base = seq[i]
        if base == 'A' or base == 'T' or  base == 'C' or base == 'G':
            pass
        else:
            x += 1
    return x



parser = argparse.ArgumentParser(description='Take the path for the gene name and input and output files')
parser.add_argument('-i', '--geneIDs', help="File containing GenBank genome IDs with one ID per line", required = True)
parser.add_argument('-g', '--geneNames', help="File containing name of gene, minsize, maxsize, and max # of Ns. Separated by whitespace and each gene on newline", required = True)

args = parser.parse_args()

geneNames_path = args.geneNames
subsamples_path = args.geneIDs
subsamples = open(subsamples_path, 'r')
db = "nucleotide"

gene_data = []
gene_handles = {}
with open(geneNames_path, 'r') as file:
    for line in file:
        geneName, minSize, maxSize, Ns = line.split()
        gene_file_name = f'{geneName}_sequences.fna'
        geneHandle = open(gene_file_name, 'w')
        gene_handles[geneName] = geneHandle
        gene_data.append({'gene_name': geneName, 'minSize': int(minSize), 'maxSize': int(maxSize), 'numN': int(Ns)})


Entrez.email = 'jd2669@msstate.edu'



for line in subsamples:
    if line[0] != '#':
        accession_number = line.strip()  # Remove newline characters
        print("Working on", accession_number)

        # Retrieve the sequence in GenBank format
        handle = Entrez.efetch(db=db, id=accession_number, rettype="gb", retmode="text")
        seq_record = SeqIO.read(handle, "genbank")
        handle.close()

        for dict in gene_data:
            geneName = dict['gene_name']
            minSize = dict['minSize']
            maxSize = dict['maxSize']
            numNs = dict['numN']
            num_printed = 0

            for feature in seq_record.features:
                keepSeq = True
                if "CDS" == feature.type and "gene" in feature.qualifiers.keys() and feature.qualifiers['gene'][0] == geneName:
                    # Create a SeqRecord for the gene sequence
                    gene_seq_record = SeqRecord(seq=feature.location.extract(seq_record.seq), id=accession_number, description="")
                    seqLen = len(gene_seq_record.seq)
                    print(f"{geneName} was found and it has a size of: {seqLen}")
                    print("Size reqs", minSize, maxSize)
                    if (minSize < 0 or  seqLen >= minSize) and (maxSize < 0 or seqLen <= maxSize):
                        if numNs > -1:
                            x = count_undef(gene_seq_record)
                            if (numNs < 1 and (x/seqLen) > numNs ) or (numNs>=1 and x > numNs) :
                                keepSeq = False

                        if keepSeq == True:
                            # Write the gene sequence to the output file
                            SeqIO.write([gene_seq_record], gene_handles[geneName], "fasta")
                            num_printed += 1
            if num_printed != 1:
                print(f'Warning, {num_printed} {geneName} sequences printed for {accession_number}')
            else:
                print(f'Successfully extracted 1 {geneName} sequence for {accession_number}')
