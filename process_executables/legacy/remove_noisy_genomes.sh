#!/bin/bash

#SBATCH --job-name=remove_noisy_genomes  ## Name of the job.
#SBATCH -A JMARTINY_LAB                ## Account to charge; personal or lab
#SBATCH -p standard                    ## Partition/queue name
#SBATCH --time=8:00:00
#SBATCH --nodes=1                      ## (-N) number of nodes the job will use
#SBATCH --mem-per-cpu=4GB
#SBATCH --cpus-per-task=1              ## number of CPUs to be used per task
#SBATCH -o %x.%A.%a.out             ## File to which STDOUT will be written, %j inserts jobid
#SBATCH -e %x.%A.%a.err             ## File to which STDERR will be written, %j inserts jobid
#SBATCH --array=1-1075%25               ## This value (1-N) depends on the number of samples


# Set folder paths
BASE=/dfs5/bio/abarrons/CA_curto/curto_core_genesDB
WD=$BASE/analysis_out/filtered_seqs
REFDB=$BASE/coregenes/hmmprofiles

# make a steps variable to send job in subatches for slurm array manager
# so we can wait some time for the first batch to do the moving of files and folders
# check it to macth with the corresponding number in the headers
steps=25


cd $REFDB
#proteins=$(for f in *.hmm; do printf '%s\n' "${f%.hmm}"; done)
FILE_NAME=$(ls *.hmm | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
protein="${FILE_NAME%.hmm}"


# Check if this is the first task
if [[ "$SLURM_ARRAY_TASK_ID" -eq 1 ]]; then
  echo -e "23G4B01\n35G2B04" > $WD/tmp.txt
fi

# Check if task ID is between 2 and n_steps, and introduce a delay
if [[ "$SLURM_ARRAY_TASK_ID" -ge 2 ]] && [[ "$SLURM_ARRAY_TASK_ID" -le $steps ]]; then
  sleep 10
fi

cd $WD

if [[ -s ${protein}.total.conc.fna ]]; then

  filterbyname.sh \
    in=${protein}.total.conc.fna \
    out=${protein}.total.conc.tmp.fna \
    names=tmp.txt \
    ow=t \
    include=f \
    fixjunk

    mv ${protein}.total.conc.tmp.fna ${protein}.total.conc.fna

  else
    echo "No files found for string: $protein"
fi

if [[ "$SLURM_ARRAY_TASK_ID" -lt $SLURM_ARRAY_TASK_MAX ]]; then
  now=$(date)
  echo "Job finished at: $now"
fi

# Check if this is the last task
if [[ "$SLURM_ARRAY_TASK_ID" -eq $SLURM_ARRAY_TASK_MAX ]]; then

  echo "Entering loop..."
  # Wait until all other tasks of this job array have finished
    while true; do
        # Count tasks with the same job array ID that are pending or running
        TASK_COUNT=$(squeue -h --name=remove_noisy_genomes -o "%A" | wc -l)
        echo "Current task count: $TASK_COUNT"

        # If count is 1 (just the last task), break out of the loop
        if [[ $TASK_COUNT -le 1 ]]; then
          echo "Task count is less than or equal to 1, breaking..."
          break
        fi
        echo "Sleeping..."
        # Sleep for a while before checking again
        sleep 20
    done

  echo "Exited loop..."
  rm $WD/tmp.txt
  now=$(date)
  echo "Job finished at: $now"

fi
