#!/bin/bash
### $ -N submit
### $ -m a
### $ -q bio

# take the aligned protein files and build individual trees and create HMM profiles
# use the concatenated tree as a "guide" tree to retain relative node structure

BASE=/bio/abchase/refDB/coregenes
BASEDIR=$BASE/individ

cd $BASEDIR

for f in *.total.aln
do
	output=${f%.total.aln}

	echo "#!/bin/bash
#$ -N ${output}.tre
#$ -m a
#$ -q bio
#$ -pe openmp 4

module load RAxML/8.2.10
module load hmmer/3.1b2 

BASE=${BASE}
BASEDIR=${BASEDIR}
reftree=\$BASE/RAxML_bipartitionsBranchLabels.coregenes

cd \$BASEDIR

rm -f RAxML_*${output}

raxmlHPC-PTHREADS-AVX -s ${f} -m PROTGAMMAWAGF -n ${output} -x 100 -# 100 -p 4321 -f a -g \$reftree -T 4

hmmbuild ${output}.hmm ${f}

	" > $output.tre.sh

	qsub $output.tre.sh
	sleep 30
done
