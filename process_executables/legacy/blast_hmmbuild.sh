#!/bin/bash

BASE=/Volumes/Bennett_BACKUP/Research/curtobacterium/curto-elevation/survey_isolates/miga_clades/roary/coregenes
OUTDIR=$BASE/referenceDB 

mkdir -p $OUTDIR

# use the final concat file to subset genomes and grab taxIDs from mappingFile
cd $BASE/master_tree

print-fasta-id.py concat.aligned.coregenes.fa

# double sure it matches up
cat concat.aligned.coregenes.fa | grep '>' | wc -l
echo "should be 75 genomes"
echo ""

rm -f $OUTDIR/allgenomeIDs.txt

# need to get NCBI taxIDs and output to mapping reference file
while read line
do

	if [[ $line == Curtobacterium_* ]]
	then
		taxid=2034

	elif [[ $line == Frigoribacterium_* ]]
	then
		taxid=96492

	else
		echo "check ${line}...something is up..."
	fi

	echo -e "${taxid}\t${line}" >> $OUTDIR/allgenomeIDs.txt

done < concat.aligned.coregenes_ids.txt


# now, if you are done with that and are satisfied, can move on and comment above out...
# save the resulting final file as $OUTDIR/allgenomeIDs.txt
# we will use this to make the BLAST DB mapping file and rename sequences to make into BLAST DBs

wc -l $OUTDIR/allgenomeIDs.txt
refgen=$OUTDIR/allgenomeIDs.txt

# make newnames for each fasta file with new headers

cd $BASE
rm -f $OUTDIR/*.blast.fasta 
rm -f $OUTDIR/FINAL_markers.map

# format for newnames == "gnl|BACT|$taxID[i]_$protein taxon=$taxID, $genomeID"


for f in *.total.fa
do
	print-fasta-id.py $f 

	protein=${f%.total.fa}

	echo "processing ${protein}..."

	count=1
	rm -f temp.txt

	while read line
	do
		refgenome=$(cat $refgen | grep "\b${line}\b")

		genomeID=$(echo ${refgenome} | cut -f2 -d' ')
		taxID=$(echo ${refgenome} | cut -f1 -d' ')

		echo "gnl|BACT|${taxID}.${count}_${protein}:taxon=${taxID},:${genomeID}" >> temp.txt
		echo "gnl|BACT|${taxID}.${count}_${protein}:${taxID}" | tr ':' '\t' >> $OUTDIR/FINAL_markers.map

		count=`expr $count + 1`

	done < $protein.total_ids.txt

	paste $protein.total_ids.txt temp.txt > newnames.$protein.txt

	# now rename the fasta files and output for BLAST DB
	fasta-rename.py $f newnames.$protein.txt $protein.blast.temp.fasta 
	cat $protein.blast.temp.fasta | tr ':' ' ' > $OUTDIR/$protein.blast.fasta 

	rm -f $protein.total_ids.txt
	rm -f temp.txt
	rm -f newnames.$protein.txt
	rm -f $protein.blast.temp.fasta

done


cd $OUTDIR

cat *.blast.fasta > total_coregenes.faa
rm *.blast.fasta 

# create a custom local database
makeblastdb -in total_coregenes.faa -dbtype 'prot' -out total_coregenes -parse_seqids -taxid_map FINAL_markers.map

# check to make sure it worked
blastdbcmd -db total_coregenes -entry all -outfmt "%T"

# should get a long list of tax ID


# now build HMMer profiles for DB as well...

# did this when I build trees on the HPC with reftrees.sh


