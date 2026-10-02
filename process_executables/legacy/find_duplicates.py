#!/usr/bin/env python3
import sys
from Bio import SeqIO
from collections import defaultdict

def find_duplicated_sequences(fasta_file):
    sequences = defaultdict(list)

    # Load sequences into the dictionary
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Convert sequence to string, use as key and append the ID as value
        sequences[str(record.seq)].append(record.id)

    # Filter out unique sequences and return only duplicates
    duplicated = {seq: ids for seq, ids in sequences.items() if len(ids) > 1}

    return duplicated

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python find_duplicates.py <path_to_fasta_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    dupes = find_duplicated_sequences(fasta_file)

    print(f"{fasta_file}|{len(dupes)}")

    for index, ids in enumerate(dupes.values(), start=1):
        print(f"Sequence {index}|{'|'.join(ids)}")
