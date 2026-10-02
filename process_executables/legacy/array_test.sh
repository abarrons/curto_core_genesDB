#!/bin/bash

#SBATCH --job-name=array_test  ## Name of the job.
#SBATCH -A JMARTINY_LAB                ## Account to charge; personal or lab
#SBATCH -p standard                    ## Partition/queue name
#SBATCH --time=36:00:00
#SBATCH --nodes=1                      ## (-N) number of nodes the job will use
#SBATCH --mem-per-cpu=4GB
#SBATCH --cpus-per-task=1              ## number of CPUs to be used per task
#SBATCH -o %x.%A.%a.out             ## File to which STDOUT will be written, %j inserts jobid
#SBATCH -e %x.%A.%a.err             ## File to which STDERR will be written, %j inserts jobid
#SBATCH --array=1-4               ## This value (1-N) depends on the number of samples

#This script will run an in-house python script to find duplicted sequences in each file and print a result of it

# Set folder paths
BASE=/dfs5/bio/abarrons/CA_curto/curto_core_genesDB
WD=$BASE/all_genomes4DB/hmmscan_out
REFDB=$BASE/coregenes/hmmprofiles
PROGRAM_PATH=$BASE/process_executables
OUTDIR=$BASE/analysis_out_rd
OUT_filter=$BASE/all_genomes4DB/filtered_seqs
# make a steps variable to send job in subatches for slurm array manager
# so we can wait some time for the first batch to do the moving of files and folders
# check it to macth with the corresponding number in the headers
steps=100


#get the list of proteins
cd $REFDB
#proteins=$(for f in *.hmm; do printf '%s\n' "${f%.hmm}"; done)
FILE_NAME=$(ls *.hmm | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
protein="${FILE_NAME%.hmm}"


cd $OUT_filter

# Check if this is the first task
if [[ "$SLURM_ARRAY_TASK_ID" -eq 1 ]]; then

  
fi


# Check if task ID is between 2 and 20, and introduce a delay
if [ "$SLURM_ARRAY_TASK_ID" -ge 2 ] && [ "$SLURM_ARRAY_TASK_ID" -le $steps ]; then
  sleep 10
fi


#load python to use pyhton specific packages
module load python/3.10.2

# Iterate over the protein names to make an file for each protein
# Find files matching the string
#files=$(find . -name "*.${protein}.*" | cut -c 3-)

# Concatenate sequences from different files into a single multi-sequence FASTA file
#cat $files > $OUTDIR/$protein.total_conc.fna


# Perform analysis with "find_duplicates.py" to generate an output that we can parse
if [[ -s $OUT_filter/${protein}.total.conc.fna ]]; then

  python $PROGRAM_PATH/legacy/find_duplicates.py $OUT_filter/${protein}.total.conc.fna > $OUTDIR/${protein}.dupes.result.txt
  # get total of duplicated sequences
  total_dupes=$(cat $OUTDIR/${protein}.dupes.result.txt | head -n 1 | cut -d "|" -f2)
  # get total of genomes with duplicated sequences for these gene
  total_genomes=$(cat $OUTDIR/${protein}.dupes.result.txt | tail -n +2 | cut -d "|" -f 2- | tr "|" "\n" | sort | uniq | wc -l)
  # get a list of sequence genomes names to filter out if necesary
  cat $OUTDIR/${protein}.dupes.result.txt | tail -n +2 | cut -d "|" -f 2- | tr "|" "\n" | sort | uniq > $OUTDIR/${protein}.dupes_seqIDs.txt
  # get count of sequences in files
  total_sequences=$(grep -c ">" $OUT_filter/${protein}.total.conc.fna)
  # get a list of the genomes where this gene was found
  grep ">" $OUT_filter/${protein}.total.conc.fna | cut -b 2- | cut -d " " -f1 > $OUTDIR/${protein}.genomes_list.txt
  # print output variables to summary file
  echo -e "$protein\t$total_sequences\t$total_dupes\t$total_genomes" > $OUTDIR/${protein}.summary_table.txt
else
  echo "No files found for string: $protein"
fi




# Check if this is the last task
if [[ "$SLURM_ARRAY_TASK_ID" -eq $SLURM_ARRAY_TASK_MAX ]]; then

  # Wait until all other tasks of this job array have finished
    while true; do
        # Count tasks with the same job array ID that are pending or running
        TASK_COUNT=$(squeue -h --name=analyze_core_genes -o "%A" | wc -l)

        # If count is 1 (just the last task), break out of the loop
        if [[ $TASK_COUNT -le 1 ]]; then
            break
        fi

        # Sleep for a while before checking again
        sleep 60
    done

  sleep 30

   # move to filtered_seqs folders
   cd $OUTDIR

   # create output file for summary results
   echo  -e "gene\ttotal_sequences\tduplicated_seqs\tn_genomes" > $OUTDIR/final_summary_table.txt
   # concatenate all summary_tables into a single fille with all genes info
   cat *.summary_table.txt >> final_summary_table.txt
   rm *.summary_table.txt
fi
