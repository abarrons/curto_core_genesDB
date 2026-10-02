#!/usr/bin/env python3
from Bio import Entrez

def get_taxid(genus):
    Entrez.email = "abarrons@uci.edu"  # Always provide your email when using NCBI's API
    handle = Entrez.esearch(db="Taxonomy", term=genus)
    record = Entrez.read(handle)
    handle.close()

    # Check if there's at least one result and return the first taxid
    if record["Count"] > "0":
        return record["IdList"][0]
    else:
        return None

genera = ["Agreia", "Bacillus", "Clavibacter", "Curtobacterium", "Frigoribacterium", "Frondihabitans", "Herbiconiux", "Leifsonia", "Leucobacter", "Microbacterium", "Mycobacterium", "Plantibacter", "Pseudoclavibacter", "Pseudomonas", "Rathayibacter", "Sphingomonas"]  # Add your list of bacterial genera here

for genus in genera:
    taxid = get_taxid(genus)
    if taxid:
        print(f"{genus}\t{taxid}")
    else:
        print(f"{genus}\tNo taxid found")
