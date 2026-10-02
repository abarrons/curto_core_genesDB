#!/usr/bin/env python3

"""
Script to calculate gene lengths for each genome from FASTA files.

This script processes FASTA files where each file corresponds to a gene and contains
the sequences for all genomes that have that gene. It calculates the length of each
gene for each genome and generates a matrix where rows are genes, columns are genomes,
and values are the lengths of the genes in the corresponding genomes (0 if the gene is absent).

Inputs:
- Directory containing FASTA files (.clean.fna format).

Outputs:
- CSV file ("gene_lengths_by_genome.csv") containing the gene lengths matrix.

Requirements:
- Python 3
- Biopython library (install with `pip install biopython`)

Usage:
    python calculate_gene_lengths.py --fasta_dir /path/to/fasta/files

Author: Alberto Barron Sandoval with chatGPT assistance
Date: 12/13/2024
"""

import os
import pandas as pd
from Bio import SeqIO
import argparse

def main(fasta_dir, output_file, count_file=None):
    # Dictionary to store gene lengths for each genome
    gene_genome_lengths = {}

    # Set to track processed genomes
    processed_genomes = set()

    # Loop through each FASTA file in the directory
    for fasta_file in os.listdir(fasta_dir):
        if fasta_file.endswith(".clean.fna"):  # Process only the relevant FASTA files
            gene_name = fasta_file.split(".")[0]  # Extract gene name
            filepath = os.path.join(fasta_dir, fasta_file)

            # Parse the FASTA file
            for record in SeqIO.parse(filepath, "fasta"):
                genome = record.id  # Extract genome name from FASTA header
                length = len(record.seq)  # Calculate sequence length

                # Track processed genomes
                processed_genomes.add(genome)

                # Add to gene-genome dictionary
                if gene_name not in gene_genome_lengths:
                    gene_genome_lengths[gene_name] = {}
                gene_genome_lengths[gene_name][genome] = length

    # Convert the dictionary into a DataFrame
    gene_length_df = pd.DataFrame.from_dict(gene_genome_lengths, orient="index").fillna(0)

    # Save to a CSV file
    gene_length_df.to_csv(output_file)
    print(f"Gene lengths by genome saved to '{output_file}'.")

    # Debugging: Compare processed genomes with count table (if provided)
    if count_file:
        counts = pd.read_csv(count_file, index_col=0)
        genomes_in_counts = set(counts.index.tolist())

        # Identify missing genomes
        missing_in_lengths = genomes_in_counts - processed_genomes
        missing_in_counts = processed_genomes - genomes_in_counts

        if missing_in_lengths:
            print("Genomes present in count table but missing in FASTA processing:")
            print(missing_in_lengths)
        if missing_in_counts:
            print("Genomes processed from FASTA but not in count table:")
            print(missing_in_counts)

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Calculate gene lengths by genome from FASTA files.")
    parser.add_argument("--fasta_dir", required=True, help="Path to directory containing FASTA files.")
    parser.add_argument("--output_file", default="gene_lengths_by_genome.csv", help="Output CSV file.")
    parser.add_argument("--count_file", help="Path to count table CSV for debugging.")
    args = parser.parse_args()

    # Run the main function
    main(args.fasta_dir, args.output_file, args.count_file)
