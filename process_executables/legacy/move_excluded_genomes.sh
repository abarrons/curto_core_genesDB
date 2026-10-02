#!/bin/bash


## Set environmental variables
BASE=/dfs5/bio/abarrons/CA_curto/curto_core_genesDB
WD=$BASE/all_genomes4DB




## move >99% ANI genomes to another folder

# make an direcotry to move excluded files
mkdir $WD/excluded_genomes
cd $WD


for file in `cat $BASE/isolates_to_remove.txt | cut -f2`; do
  echo "Moving $file"
  # conditianal to check if any files for that genomes where found
  if [[ -s $file  ]]; then
    # if files are found the move the exclude folder
    mv $file $WD/excluded_genomes

  else
    # if files are not found then just print a message
    echo "No files where found for file: $file"
  fi

done
