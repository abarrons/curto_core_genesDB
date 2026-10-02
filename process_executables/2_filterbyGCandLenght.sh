#!/bin/bash

#SBATCH --job-name=filterbylength&gc		  ## Name of the job.
#SBATCH -A JMARTINY_LAB             ## Account to charge; personal or lab
#SBATCH -p standard                 ## Partition/queue name
#SBATCH --nodes=1                   ## (-N) number of nodes the job will use
#SBATCH --ntasks=1                  ## (-n) number of processes to be launched
#SBATCH --cpus-per-task=4           ## number of OpenMP threads - total RAM request = 4 * 4.5GB/core = 18GB
#SBATCH -o filterby_%j.out          ## File to which STDOUT will be written, %j inserts jobid
#SBATCH -e filterby_%j.err          ## File to which STDERR will be written, %j inserts jobid
#SBATCH --array=1-37               ## This value (1-N) depends on the number of samples

##set environmental variables
export BASE=/dfs5/bio/abarrons/CA_curto/curto_core_genesDB
export WD=$BASE/blobtools_out

#Get sequence IDs from blast results that will be filtered for the next analysis
## go to blast results directory
cd $WD
ISOLATE_ID=`ls *.table.txt | cut -d "_" -f 1 | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
mkdir seqs_ids

tail -n +12 ${ISOLATE_ID}_*.table.txt | awk ' $2 > 5000 && $5 > 30' | cut -f1 > seqs_ids/${ISOLATE_ID}.tmp.txt

#load BBmap suite
module load bbmap/38.96

cd $BASE/og_genomes

# filter each genome using BBMap software
filterbyname.sh \
 in=${ISOLATE_ID}_scaffolds.fasta \
 out=${ISOLATE_ID}_filtered_scaffolds.fasta \
 names=$WD/seqs_ids/${ISOLATE_ID}.tmp.txt \
 ow=t \
 include=t \
 fixjunk

# remove unflitered assemblies
rm ${ISOLATE_ID}_scaffolds.fasta

## rename file with isolate id
# loop to rename files
cat  $BASE/og_isolates_b* | grep ${ISOLATE_ID} | while IFS=$'\t' read -r isID seqID tax_id blah file; do
  # Remove any leading/trailing whitespaces from the new_name
  new_name=$(echo "$isID" | tr -d '[:space:]')
  # Print a debug message with old and new names
  echo "Renaming '${seqID}_filtered_scaffolds.fasta' to '${new_name}.fna'"
    # copy the file using 'cp' command if it exists in the directory
    if [ -e "${seqID}_filtered_scaffolds.fasta" ]; then
         mv "${seqID}_filtered_scaffolds.fasta" "${new_name}.fna"
    else
        # Print a warning if the file is not found
        echo "Warning: File '${seqID}_filtered_scaffolds.fasta' not found. Skipping renaming."
    fi
done
