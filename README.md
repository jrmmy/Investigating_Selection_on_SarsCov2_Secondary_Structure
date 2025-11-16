# Pipeline README

This repository contains scripts for downloading, filtering, aligning, trimming, and analyzing gene sequences using a combination of **NCBI Entrez**, **Biopython**, **TranslatorX**, **IQ-TREE**, **HyPhy**, and SLURM batch jobs.

It consists of:

* A **SLURM job** that builds per-gene analysis pipelines.
* A **frame validator** that removes indel-disrupted ORFs.
* A **GenBank gene extractor** that downloads genomes and extracts specified genes.

Below is a detailed explanation of each script, inputs, outputs, and workflow.

---

## 1. `alignAll.sh`

A Bash script that:

* Reads a **guidelines file** defining gene_name, min_size, max_size, and maximum Ns.
* For each gene:

  * Creates an output folder.
  * Creates a **SLURM job script** that performs:

    * Amino-acid guided nucleotide alignment using **TranslatorX**.
    * Phylogenetic tree construction with **IQ-TREE**.
    * Creation of a PAML control file (optional if using codeml).
    * Alignment trimming via custom Python script (`trim.py`).
    * Selection analysis using **HyPhy SLAC**.
    * Conversion of SLAC output to CSV.

### Usage

```
./run_guidelines.sh guidelines.txt input_sequences_folder output_folder baseline_codeml.ctl
```

### Inputs

* `guidelines.txt` — whitespace‑separated list:

  ```
  geneName   minSize   maxSize   maxN
  ```
* `input_sequences_folder` — contains `<geneName>_sequences.fna` files.
* `output_folder` — where results will be stored.
* `baseline_codeml.ctl` — template codeml file with placeholders.

### Outputs (per gene)

Inside `<geneName>_output/`:

* `<geneName>.nt_ali.fasta` — aligned nucleotide sequences.
* `<geneName>.treefile` — IQ‑TREE phylogeny.
* `<geneName>.fna.trimmed` — trimmed alignment.
* `<geneName>.slac` — HyPhy SLAC results.
* `<geneName>.slac.csv` — parsed CSV output.

---

## 2. `FindDashes.py`

Ensures that aligned genomes remain **codon-aligned** (no frameshift‑causing indels). It removes sequences where gaps occur in non‑triplet multiples.

### How it works

* Parses a FASTA file.
* Checks each sequence for gap runs of length % 3 == 0.
* Writes only valid sequences to the output FASTA.

### Usage

```
python3 frame_validator.py -i aligned_input.fasta -o valid_output.fasta
```

### Output

* A FASTA file containing only codon-preserving sequences.
* Terminal summary of valid vs rejected sequences.

---

## 3. `ChopByGene.py`

Downloads GenBank records from NCBI and extracts CDS features for selected genes.

### What it does

* Reads:

  * A file of **accession numbers** (one per line).
  * A **gene config file** listing geneName, minSize, maxSize, maxNs.
* For each genome:

  * Fetches the GenBank file via **Entrez.efetch**.
  * Extracts CDS features whose gene qualifier matches the target gene.
  * Applies filtering:

    * Sequence length must be within minSize / maxSize (if positive).
    * Number of ambiguous bases must be ≤ maxNs.
  * Writes successful gene sequences to: `<geneName>_sequences.fna`.

### Usage

```
python extract_genes.py -i genome_IDs.txt -g gene_config.txt
```

### Example gene_config.txt

```
Nsp1   250   300   0
S      3500  4000  5
N      -1    -1    -1
```

### Outputs

* For each gene, a file named:

  * `<geneName>_sequences.fna`
* Log messages reporting:

  * Whether gene was found.
  * Sequence size.
  * Whether it met filtering criteria.
  * If more or fewer than 1 sequence was extracted.

---

## 4. Dependencies

This workflow requires the following software:

### Modules / External Tools

* TranslatorX
* Seaview
* IQ-TREE (v2.0 or later)
* HyPhy
* (Optional) PAML

### Python

* Biopython
* argparse
* csv

### Conda environments

Used in the SLURM script:

* A Biopython environment for trimming.
* A HyPhy environment for SLAC execution.

---

## 5. Workflow Overview

1. **Extract target genes** from full genomes using `extract_genes.py`.
2. **Validate ORFs** in the aligned sequences using `frame_validator.py`.
3. **Prepare SLURM jobs** using `run_guidelines.sh`.
4. For each gene, SLURM will:

   * Align sequences.
   * Build a phylogeny.
   * Trim alignments.
   * Run HyPhy SLAC.
   * Convert results to CSV.

---

## 6. Notes

* Ensure your Entrez email is updated inside the gene extraction script.
* If using on a cluster, uncomment the `sbatch` command to auto‑submit jobs.
* Codeml support is present but commented out; enable if needed.

