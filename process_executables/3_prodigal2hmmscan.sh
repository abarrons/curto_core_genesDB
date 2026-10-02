#!/bin/bash

#SBATCH --job-name=prodigal2hmmscan  ## Name of the job.
#SBATCH -A JMARTINY_LAB                ## Account to charge; personal or lab
#SBATCH -p standard                    ## Partition/queue name
#SBATCH --time=12:00:00
#SBATCH --nodes=1                      ## (-N) number of nodes the job will use
#SBATCH --cpus-per-task=8               ## number of CPUs to be used per task
#SBATCH --mem-per-cpu=4GB
#SBATCH -o %x.%A.%a.out             ## File to which STDOUT will be written, %j inserts jobid
#SBATCH -e %x.%A.%a.err             ## File to which STDERR will be written, %j inserts jobid
#SBATCH --array=1-470                   ## This value (1-64) depends on the number of samples (1-no. of samples)

## This sets environmnetal variables to be used by slurm in order to use all possible RAM
## Set environmental variables
BASE=/dfs5/bio/abarrons/CA_curto/curto_core_genesDB
REFDB=$BASE/coregenes/hmmprofiles
WD=$BASE/all_genomes4DB
OUT_prodigal=$WD/translated_seqs
OUT_hmmr=$WD/hmmscan_out
OUT_filter=$WD/filtered_seqs
EXCLUDE_LIST=/dfs5/bio/abarrons/CA_curto/isolates/all_good_genomes/checkM_out/plus5pcnt_contam_genomes.txt
THREADS=8

# Check if this is the first task
if [[ "$SLURM_ARRAY_TASK_ID" -eq 1 ]]; then

  # make a clean up of the folder and create output directories
  rm -rf $OUT_prodigal
  mkdir $OUT_prodigal
  rm -rf $OUT_hmmr
  mkdir $OUT_hmmr
  mkdir $OUT_hmmr/seqs_ids
  rm -rf $OUT_filter
  mkdir $OUT_filter

  # move >99% ANI genomes to another folder

  # make an direcotry to move excluded files
  #rm -rf $WD/excluded_genomes
  #mkdir $WD/excluded_genomes

  #for file in `cat $BASE/isolates_to_remove.txt | cut -f2`; do
  #  echo "Moving $file"
  #  # conditianal to check if any files for that genomes where found
  #  if [[ -s $file  ]]; then
  #    # if files are found the move the exclude folder
  #    mv $file $WD/excluded_genomes

  #  else
  #    # if files are not found then just print a message
  #    echo "No files where found for file: $file"
  #  fi

  #done

  cd $WD

  for pattern in $(cat $EXCLUDE_LIST); do
    echo "Searching for files matching pattern: $pattern"

    # Find files that match the pattern
    matching_files=$(find . -type f -name "*${pattern}*")

    # Check if any files were found
    if [[ -n "$matching_files" ]]; then
        echo "Moving files matching pattern: $pattern"
        # Move each matching file
        for file in $matching_files; do
            mv "$file" $WD/excluded_genomes
        done
    else
        # If no files match the pattern, print a message
        echo "No files found for pattern: $pattern"
    fi
done

fi

# Check if task ID is between 2 and 20, and introduce a delay
if [[ "$SLURM_ARRAY_TASK_ID" -ge 2 ]]; then
  sleep 60
fi


## Load programs: prodigal
module load prodigal/2.6.3

## Go to working directory and make folder for output files
cd $WD


## Set environmental variables specific for the job
SAMPLE_NAME=`ls *.fna | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
ISOLATE_ID=`basename $SAMPLE_NAME .fna`


## Prodigal will translate nucleotide sequences to protein sequences for subsequent analysis
prodigal \
  -a $OUT_prodigal/${ISOLATE_ID}.faa \
  -d $OUT_prodigal/${ISOLATE_ID}_nucleotid_seqs.fna \
  -i $SAMPLE_NAME


## Load programs: HMMER
module load hmmer/3.3

#get the list of proteins
cd $REFDB
proteins=$(for f in *.hmm; do printf '%s\n' "${f%.hmm}"; done)

## Go to working directory and make folder for output files
cd $OUT_prodigal


### hmmscan will identify the protien sequences that mache to the core gene sequences in the .hmm profiles (downloaded from Alex's github: elevation_community)
### make a while read line loop to go through each gene

for protein in $proteins
  do
        hmmsearch --tblout $OUT_hmmr/${ISOLATE_ID}.${protein}.hmm.txt -E 1e-10 --cpu $THREADS \
              $REFDB/${protein}.hmm ${ISOLATE_ID}.faa
  done


#load BBmap suit we eill need it to filter the fasta files
module load bbmap/38.96


## Go to working directory and make folder for output files
cd $OUT_hmmr


## make a loop to go throughg all the combination of genomes and proteins
## to extract the sequence identifiers for each protein in each genome and
## get that sequence from the protein fasta file of that genome

#first go to the hmmscan outp folder

for protein in $proteins
  do

    # first get the highest socore match for each protein and store the protein name and sequenceID
    # to filter that sequence from the nucleotide fasta file
      hmmscanout=$(cat ${ISOLATE_ID}.${protein}.hmm.txt | awk '$1 !~ /^#/ && ($6 > max_score || NR == 1) { max_score = $6; best_entry = $0 } END { split(best_entry, fields); prot_name = fields[3]; seqID = fields[1]; print prot_name, seqID }')

    # Extract the variables from the output
      prot_name=$(echo "$hmmscanout" | awk '{print $1}')
      seqID=$(echo "$hmmscanout" | awk '{print $2}')

      # Check if $seqID is not empty
      if [[ -z "$seqID" ]]; then
        # print warning message that this gene was not foun in this geneome
        echo "Warning: gene $protein not found in genome: $ISOLATE_ID"
        # store that info in a summary file
        echo  -e "${ISOLATE_ID}\t${protein}" >> $OUT_filter/${ISOLATE_ID}.missing_genes.txt
        continue
      fi

      echo "$seqID" > $OUT_hmmr/seqs_ids/${ISOLATE_ID}.${protein}.tmp.txt


    # now let's use those variables to get the mathcing sequence from each genome
    # we will use filterbyname fromm the BBmap tool suite
      filterbyname.sh \
        in=$OUT_prodigal/${ISOLATE_ID}_nucleotid_seqs.fna \
        out=$OUT_filter/${ISOLATE_ID}.${prot_name}.fna \
        names=$OUT_hmmr/seqs_ids/${ISOLATE_ID}.${protein}.tmp.txt \
        ow=t \
        include=t \
        fixjunk


    # Replace header in the FASTA file
      fasta_file="$OUT_filter/${ISOLATE_ID}.${prot_name}.fna"
      tmp_file="$OUT_filter/${ISOLATE_ID}.${prot_name}.temp.fna"

   # Use awk to replace specific parts of the header
   awk -v var1="$ISOLATE_ID" '/^>/{ sub(">[^ ]*", ">" var1); } 1' "$fasta_file" > "$tmp_file"

   # Rename the temporary file to overwrite the original FASTA file
   mv "$tmp_file" "$fasta_file"


 done


# Check if this is the last task
if [[ "$SLURM_ARRAY_TASK_ID" -eq $SLURM_ARRAY_TASK_MAX ]]; then

  echo "Entering loop..."
  # Wait until all other tasks of this job array have finished
    while true; do
        # Count tasks with the same job array ID that are pending or running
        TASK_COUNT=$(squeue -h --name=prodigal2hmmscan -o "%A" | wc -l)
        echo "Current task count: $TASK_COUNT"

        # If count is 1 (just the last task), break out of the loop
        if [[ $TASK_COUNT -le 1 ]]; then
          echo "Task count is less than or equal to 1, breaking..."
          break
        fi

        # Sleep for a while before checking again
        echo "Sleeping..."
        sleep 60
    done

  echo "Exited loop..."
  sleep 30

  rm -rf $OUT_hmmr



   # move to filtered_seqs folders
   cd $OUT_filter

   for protein in $proteins
     do

       tmp_file=$(mktemp)
       cat *.${protein}.total.fna > $tmp_file

       if [[ -s $tmp_file ]]; then

         rm *.${protein}.total.fna
         echo "concatenating files for gene $protein"
         cat $tmp_file > ${protein}.total.conc.fna
         rm $tmp_file

       else
         echo "No files found for gene $protein"
       fi
     done

     # now concatenate all missing gene info into a single file
     # create a file to store info about which genes were not found and the genome where it was not found.
     echo "Concatenating missing gene summary tables..."
     echo -e "genome\tgene" > $OUT_filter/summary_missing_genes.txt
     cat *.missing_genes.txt >> summary_missing_genes.txt
     rm *.missing_genes.txt

   sleep 60
   echo "Sending next job: 4_derep_seqs.sh"
   cd $BASE/process_logs
   sbatch $BASE/process_executables/4_derep_seqs.sh
   now=$(date)
   echo "Job finished at: $now"
fi
