#!/bin/bash

#SBATCH --job-name=move_og_files  ## Name of the job.
#SBATCH -A JMARTINY_LAB                ## Account to charge; personal or lab
#SBATCH -p standard                    ## Partition/queue name
#SBATCH --time=24:00:00
#SBATCH --nodes=1                      ## (-N) number of nodes the job will use
#SBATCH --mem=4GB                     ## number of OpenMP threads - total RAM request = 4 * 4.5GB/core = 18G
#SBATCH --cpus-per-task=1              ## number of CPUs to be used per task
#SBATCH -o move_og_files_%j.out             ## File to which STDOUT will be written, %j inserts jobid
#SBATCH -e move_og_files_%j.err             ## File to which STDERR will be written, %j inserts jobid

# This script will rename the filtered_scaffolds to remove the sequence ID i.e. P002 to the actual ISOLATE_ID
# For this we need to prepare a mapping_file that has the seqIDs and the matching ISOLATE_IDs in a single line


# set environmental variables for the base and working directories
BASE=/dfs5/bio/abarrons/CA_curto/isolates
DB=/dfs5/bio/abarrons/CA_curto/curto_core_genesDB
B1=$BASE/spades_out/assemblies
BP1=$BASE/blobtools_out2
B2=$BASE/2ndbatch/spades_out/assemblies
BP2=$BASE/2ndbatch/blobtools_out

# Define the mapping file path
mapping_file1="$DB/og_isolates_b1_md.txt"
mapping_file2="$DB/og_isolates_b2_md.txt"
# Go to the directory containing the files with old names
cd $DB
# make directoryfor blobplot files and genomes
mkdir og_genomes
mkdir blobtools_out


# loop to move files from batch 1
cat  "$mapping_file1" | while IFS=$'\t' read -r isID seqID tax_id blah file; do
  # Print a debug message with old and new names
  echo "Moving '${file}' to '${DB/og_genomes}'"
    # copy the file using 'cp' command if it exists in the directory
    if [ -e "$B1/${file}" ]; then
         cp $B1/${file} $DB/og_genomes
    else
        # Print a warning if the file is not found
        echo "Warning: File '${file}' not found. Skipping renaming."
    fi
    # now lets move the blobplot files
    cp $BP1/${seqID}* $DB/blobtools_out
done


# loop to move files from batch 2
cat  "$mapping_file2" | while IFS=$'\t' read -r isID seqID tax_id blah file; do
  # Print a debug message with old and new names
  echo "Moving '${file}' to '${DB/og_genomes}'"
    # copy the file using 'cp' command if it exists in the directory
    if [ -e "$B2/${file}" ]; then
         cp $B2/${file} $DB/og_genomes
    else
        # Print a warning if the file is not found
        echo "Warning: File '${file}' not found. Skipping renaming."
    fi
    # now lets move the blobplot files
    cp $BP2/${seqID}* $DB/blobtools_out
done
