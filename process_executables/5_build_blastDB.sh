#!/bin/bash

#SBATCH --job-name=buildDB  ## Name of the job.
#SBATCH -A JMARTINY_LAB                ## Account to charge; personal or lab
#SBATCH -p standard                    ## Partition/queue name
#SBATCH --time=2:00:00
#SBATCH --nodes=1                      ## (-N) number of nodes the job will use
#SBATCH --mem=4GB                     ## number of OpenMP threads - total RAM request = 4 * 4.5GB/core = 18G
#SBATCH --cpus-per-task=1              ## number of CPUs to be used per task
#SBATCH -o %x.%A.%a.out             ## File to which STDOUT will be written, %j inserts jobid
#SBATCH -e %x.%A.%a.err             ## File to which STDERR will be written, %j inserts jobid


# Set folder paths
BASE=/dfs5/bio/abarrons/CA_curto/curto_core_genesDB
REFDB=$BASE/coregenes/hmmprofiles
WD=$BASE/clean_genes
PROGRAM_PATH=/dfs5/bio/abarrons/programs
OUTDIR=$BASE/referenceDB
#GENES4DB=$BASE/DBgenes.txt

## Set environmental variables specific for the job
#gene=`cat $GENES4DB | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

#get the list of proteins
cd $REFDB
GENES4DB=$(for f in *.hmm; do printf '%s\n' "${f%.hmm}"; done)



rm -rf $OUTDIR
mkdir -p $OUTDIR

# verify that files listed in the the masterMD4DB.txt file are exactly the same contained in the genomes folder
# get count of files listed in masterMD4DB.txt
echo "We are expecting 470 files"
echo ""
filesinfile=$(cut -f2 $BASE/masterMD4DB.txt | sort | uniq | wc -l)
echo "##########################"
echo ""
echo "# of files listed in masterMD4DB.txt: $filesinfile"
filesindir=$(ls $BASE/all_genomes4DB/*.fna | sort | uniq | wc -l)
echo "##########################"
echo ""
echo "# of files listed in all_genomes4DB folder: $filesindir"
echo ""

# next let's compare both lists to check that all files match otherwise
# finish the script and print an error message so we can check for error or missmatches
listcomm=$(comm -3 <(ls $BASE/all_genomes4DB/*.fna | cut -d "/" -f8 | sort) <(cut -f2 $BASE/masterMD4DB.txt | sort) | wc -l)

# use a conditional to check that both lists are exactly the same
# otherwise finish the job and print error message
if [[ $listcomm -ne 0 ]]; then
  echo "Files don't match! Verify that files match and retry. Exiting..."
  exit 1
else
  echo "Files match, proceedind to build DB"
fi


rm -f $OUTDIR/allgenomeIDs.txt
# get NCBI taxids from a table that matches taxon to file id
taxidref=$BASE/ncbi_taxids.txt

# need to get NCBI taxIDs and output to mapping reference file
while read taxon; do

  reftaxon=$(cat $taxidref | grep "\b${taxon}\b")
  taxid=$(echo $reftaxon | cut -d " " -f2)

  echo "$taxid" >> $OUTDIR/allgenomestaxids.text

done < <(cut -f1 $BASE/masterMD4DB.txt)

paste $BASE/masterMD4DB.txt $OUTDIR/allgenomestaxids.text > $OUTDIR/allgenomeIDs.txt

rm -f $OUTDIR/allgenomestaxids.text

# now that we added the taxid variable to the mapping file we can proceed
# with renaming the sequence IDs in each protein fasta file

wc -l $OUTDIR/allgenomeIDs.txt
refgen=$OUTDIR/allgenomeIDs.txt

## make newnames for each fasta file with new headers
# move to folder with derep seqs
cd $WD

# format for newnames == "gnl|BACT|$taxID[i]_$protein taxon=$taxID, $genomeID"

while read gene; do

  protein=${gene}

  echo "processing ${protein}"

  #count=1

  while read genomeID; do
    echo "Processing header $genomeID for gene $protein"

    refline=$(grep -F "$genomeID" $refgen)
    taxID=$(echo $refline | cut -d " " -f6)
    correctgenomeID=$(echo $refline | cut -d " " -f5)

    if grep -q "gnl|BACT|${protein}-${correctgenomeID}:" newnames.${protein}.txt; then
      echo "Error: Duplicate genomeID ${correctgenomeID} found in newnames.${protein}.txt, printing colected variables.."
      echo "$refline"
      echo "$genomeID"
      echo "$taxID"
      echo  "$correctgenomeID"
      exit 1
    fi

    echo -e "${genomeID}\tgnl|BACT|${protein}-${correctgenomeID}:taxon=${taxID}" >> newnames.${protein}.txt
    echo "gnl|BACT|${protein}-${correctgenomeID}:${taxID}" | tr ':' '\t' >> $OUTDIR/FINAL_markers.map

    #count=`expr $count + 1`

  done < <(grep ">" ${protein}.clean.fna | cut -c 2-)

  # now rename the fasta files and output for BLAST DB
  $PROGRAM_PATH/seqkit replace -p '^(\S+)' -r '{kv}' -k newnames.${protein}.txt ${protein}.clean.fna > ${protein}.blast.temp.fasta
  cat ${protein}.blast.temp.fasta | tr ':' ' ' > $OUTDIR/${protein}.blast.fasta

  #do a clean up of intermediate files
  rm -f newnames.${protein}.txt
	rm -f ${protein}.blast.temp.fasta

done <<< $GENES4DB

module load ncbi-blast/2.13.0

cd $OUTDIR

cat *.blast.fasta > total_coregenes.fna
rm *.blast.fasta

# Check for duplicate headers
if [[ $(grep "^>" total_coregenes.fna | sort | uniq -d | wc -l) -gt 0 ]]; then
    echo "Error: Duplicate headers found in total_coregenes.fna:"
    grep "^>" total_coregenes.fna | sort | uniq -d
    exit 1
fi


# create a custom local database
makeblastdb -in total_coregenes.fna -dbtype 'nucl' -out total_coregenes -parse_seqids -taxid_map FINAL_markers.map

# check to make sure it worked
blastdbcmd -db total_coregenes -entry all -outfmt "%T"

# should get a long list of tax ID
