#!/bin/bash

#SBATCH --job-name=derep_seqs  ## Name of the job.
#SBATCH -A JMARTINY_LAB                ## Account to charge; personal or lab
#SBATCH -p standard                    ## Partition/queue name
#SBATCH --time=8:00:00
#SBATCH --nodes=1                      ## (-N) number of nodes the job will use
#SBATCH --mem-per-cpu=4GB
#SBATCH --cpus-per-task=4              ## number of CPUs to be used per task
#SBATCH -o %x.%A.%a.out             ## File to which STDOUT will be written, %j inserts jobid
#SBATCH -e %x.%A.%a.err             ## File to which STDERR will be written, %j inserts jobid
#SBATCH --array=1-1075%200               ## This value (1-N) depends on the number of samples

#This script dereplicates the sequences of each core gene with CD-HIT and summarizes prevalence stats

# Set folder paths
BASE=/dfs5/bio/abarrons/CA_curto/curto_core_genesDB
WD=$BASE/all_genomes4DB/filtered_seqs
REFDB=$BASE/coregenes/hmmprofiles
PROGRAM_PATH=/dfs5/bio/abarrons/programs/cdhit
OUTDIR=$BASE/clean_genes
THREADS=4
# make a steps variable to send job in subatches for slurm array manager
# so we can wait some time for the first batch to do the moving of files and folders
# check it to macth with the corresponding number in the headers
steps=200

initial_start=$(date +%s)

#get the list of proteins
cd $REFDB
#proteins=$(for f in *.hmm; do printf '%s\n' "${f%.hmm}"; done)
FILE_NAME=$(ls *.hmm | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
protein="${FILE_NAME%.hmm}"


cd $WD

# Check if this is the first task
if [[ "$SLURM_ARRAY_TASK_ID" -eq 1 ]]; then

  # make output directory for alignmets
  rm -rf $OUTDIR
  mkdir $OUTDIR

fi


# Check if task ID is between 2 and 20, and introduce a delay
if [[ "$SLURM_ARRAY_TASK_ID" -ge 2 ]] && [[ "$SLURM_ARRAY_TASK_ID" -le $steps ]]; then
  sleep 10
fi


# Dereplicate this gene's sequences with CD-HIT
if [[ -s ${protein}.total.conc.fna ]]; then

    ## Replace header in the FASTA file
    fasta_file="${protein}.total.conc.fna"
    tmp_file="${protein}.total.conc.tmp.fna"

    ## Use awk to remove extra info from the header
    awk '/^>/ {print $1} !/^>/ {print}' "$fasta_file" > "$tmp_file"

    ## Rename the temporary file to overwrite the original FASTA file
    mv "$tmp_file" "$fasta_file"

    # Identify and cluster duplicated sequences to have unique sequences in the
    # input fasta file with all sequnces for each gene, in a single file each
    # we will use CD-HIT for this
    $PROGRAM_PATH/cd-hit \
      -i $fasta_file \
      -o $OUTDIR/${protein}.clean.fna \
      -c 1 -g 1 -d 0 -T $THREADS -G 0 -aL 0.75 -aS 1

    # get a list of the genomes removed from this gene
    comm -23 <(grep "^>" $fasta_file | sed 's/^>//g' | sort) <(grep "^>" $OUTDIR/${protein}.clean.fna | sed 's/^>//g' | sort) > $OUTDIR/${protein}_dupes_removed_seqIDs.txt

    # get a list of curto duplicated seqs removed from this gene
    comm -12 <(cat $BASE/curtos4DBfiles.txt | rev | cut -d "." -f 2- | rev | sort) <(cat $OUTDIR/${protein}_dupes_removed_seqIDs.txt | sort) > $OUTDIR/${protein}_curto_dupes_removed_seqIDs.txt

    # get a list of curto genomes present in dereplicated file
    comm -12 <(cat $BASE/curtos4DBfiles.txt | rev | cut -d "." -f 2- | rev | sort) <(grep "^>" $fasta_file | sed 's/^>//g' | sort) > $OUTDIR/${protein}_curtosB4_derep.txt

    # get a list of Curto genomes IDs where this gene was found after dereplication
    comm -12 <(cat $BASE/curtos4DBfiles.txt | rev | cut -d "." -f 2- | rev | sort) <(grep "^>" $OUTDIR/${protein}.clean.fna | sed 's/^>//g' | sort) > $OUTDIR/${protein}_curtos_after_derep.txt

    # get count of total curto genomes that contain this gene
    total_curto=$(cat $OUTDIR/${protein}_curtosB4_derep.txt | wc -l)

    # get count of unique curto sequences for these gene after dereplication
    unique_curto_seqs=$(cat $OUTDIR/${protein}_curtos_after_derep.txt | wc -l)

    # get count of sequences in the original file, a count of how many genomes have this gene
    total_seqs=$(grep -c ">" $fasta_file)

    # get count of sequences in the dereplicated file, unique sequences
    unique_seqs=$(grep -c ">" $OUTDIR/${protein}.clean.fna)

    # calculate the percentage of unique curtos where this gene was found
    pcnt_curto_unique=$(printf "%.2f" $(echo "scale=4; ($unique_curto_seqs / 395) * 100" | bc))

    # calculate the percentage of curtos where this gene was found
    pcnt_curto_total=$(printf "%.2f" $(echo "scale=4; ($total_curto / 395) * 100" | bc))

    # calculate the percentage of unique genomes where this gene was found
    pcnt_unique=$(printf "%.2f" $(echo "scale=4; ($unique_seqs / 468) * 100" | bc))

    # calculate the percentage of all genomes where this gene was found
    pcnt_total=$(printf "%.2f" $(echo "scale=4; ($total_seqs / 468) * 100" | bc))

    # get count of removed curtos
    removed_curtos=$(cat $OUTDIR/${protein}_curto_dupes_removed_seqIDs.txt | wc -l)

    #get count of total seqs removed
    total_seqs_removed=$((total_seqs - unique_seqs))

    # print output variables to summary file
    echo -e "$protein\t$total_curto\t$unique_curto_seqs\t$unique_seqs\t$total_seqs\t$pcnt_curto_unique\t$pcnt_curto_total\t$pcnt_unique\t$pcnt_total\t$removed_curtos\t$total_seqs_removed" > $OUTDIR/${protein}.summary_table.txt
  else
    echo "No files found for string: $protein"
fi

if [[ "$SLURM_ARRAY_TASK_ID" -lt $SLURM_ARRAY_TASK_MAX ]]; then
  now=$(date)
  echo "Job finished at: $now"
  echo ""
  total_end=$(date +%s)
  total_runtime=$(echo "$total_end - $initial_start" | bc -l)
  echo "################################################################################################"
  echo ""
  echo "Total time: $total_runtime seconds"

fi


# Check if this is the last task
if [[ "$SLURM_ARRAY_TASK_ID" -eq $SLURM_ARRAY_TASK_MAX ]]; then

  echo "Entering loop..."
  # Wait until all other tasks of this job array have finished
    while true; do
        # Count tasks with the same job array ID that are pending or running
        TASK_COUNT=$(squeue -h --name=derep_seqs -o "%A" | wc -l)
        echo "Current task count: $TASK_COUNT"

        # If count is 1 (just the last task), break out of the loop
        if [[ $TASK_COUNT -le 1 ]]; then
          echo "Task count is less than or equal to 1, breaking..."
          break
        fi
        echo "Sleeping..."
        # Sleep for a while before checking again
        sleep 15
    done

  echo "Exited loop..."
  sleep 10

   # move to filtered_seqs folders
   cd $OUTDIR
   echo "Concatenating final summary table.."

   # create output file for summary results
   # concatenate all summary_tables into a single fille with all genes info
   cat *.summary_table.txt > final_summary_table.tmp.txt
   rm *.summary_table.txt
   # add header to the output file
   echo -e "protein\ttotal_curto\tunique_curto_seqs\tunique_seqs\ttotal_seqs\tpcnt_curto_unique\tpcnt_curto_total\tpcnt_unique\tpcnt_total\tn_removed_curtos\ttotal_seqs_removed" > $OUTDIR/final_summary_table.txt

   cat final_summary_table.tmp.txt | sort -t$'\t' -k3,3nr >> $OUTDIR/final_summary_table.txt
   # rmeove unsorted file
   rm final_summary_table.tmp.txt

   # print end of process messages
   now=$(date)
   echo "Job finished at: $now"
   echo ""
   total_end=$(date +%s)
   total_runtime=$(echo "$total_end - $initial_start" | bc -l)
   echo "################################################################################################"
   echo ""
   echo "Total time: $total_runtime seconds"

fi
