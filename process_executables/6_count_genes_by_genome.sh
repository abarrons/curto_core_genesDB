#!/bin/bash

## Set environmental variables
BASE=/dfs5/bio/abarrons/CA_curto/curto_core_genesDB
WD=$BASE/clean_genes
REFDB=$BASE/coregenes/hmmprofiles



#get the list of proteins
cd $REFDB
for f in *.hmm; do printf '%s\n' "${f%.hmm}"; done > $BASE/DBgenes.txt

DB_GENES_FILE=$BASE/DBgenes.txt

# Enable strict error handling
set -e
# count genes present in each genome and put it in a list


# Temporary file to store genome-to-gene data
TEMP_FILE="genome_to_genes_temp.txt"
OUTPUT_FILE="genome_gene_counts.csv"

# Move to working directory
cd "$WD"

# Ensure working directory contains files
if ! ls *_curtos_after_derep.txt 1> /dev/null 2>&1; then
  echo "No input files found in $WD"
  exit 1
fi

# Ensure DB genes file exists
if [ ! -f "$DB_GENES_FILE" ]; then
  echo "DBgenes.txt not found at $DB_GENES_FILE"
  exit 1
fi

# Read database genes into an array
declare -A db_genes
while IFS= read -r gene; do
  db_genes["$gene"]=1
done < "$DB_GENES_FILE"

# Clear temporary file if it exists
> $TEMP_FILE


for file in "$WD"/*_curtos_after_derep.txt; do
  gene_name=$(basename "$file" "_curtos_after_derep.txt")  # Extract gene name

  # Skip processing if the gene is not in the database
  if [[ -z "${db_genes[$gene_name]}" ]]; then
    continue
  fi

  # Print progress message
  echo "Processing file: $file"

    # Read each genome in the file and map to the gene
    while IFS= read -r genome; do
        echo -e "$genome\t$gene_name" >> $TEMP_FILE
    done < "$file"
done

# Sort and group by genome to count unique genes
awk -F'\t' '{
    genome_to_genes[$1][$2]=1
}
END {
    print "Genome,Gene_Count" > "'$OUTPUT_FILE'"
    for (genome in genome_to_genes) {
        print genome "," length(genome_to_genes[genome]) >> "'$OUTPUT_FILE'"
    }
}' $TEMP_FILE

# Clean up
rm $TEMP_FILE

echo "Gene counts saved to $OUTPUT_FILE"
