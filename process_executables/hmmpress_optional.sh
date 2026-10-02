#!/bin/bash

#SBATCH --job-name=hmmpress  ## Name of the job.
#SBATCH -A JMARTINY_LAB                ## Account to charge; personal or lab
#SBATCH -p standard                    ## Partition/queue name
#SBATCH --time=24:00:00
#SBATCH --nodes=1                      ## (-N) number of nodes the job will use
#SBATCH --mem=18GB                     ## number of OpenMP threads - total RAM request = 4 * 4.5GB/core = 18G
#SBATCH --cpus-per-task=1              ## number of CPUs to be used per task
#SBATCH -o hmmpress_%j.out             ## File to which STDOUT will be written, %j inserts jobid
#SBATCH -e hmmpress_%j.err             ## File to which STDERR will be written, %j inserts jobid



## Set environmental variables

REFDB=/dfs5/bio/abarrons/CA_curto/curto_core_genesDB/coregenes/hmmprofiles
## Load programs: HMMER
module load hmmer/3.3

cd $REFDB

for f in *.hmm
  do
    hmmpress -f $f
  done
