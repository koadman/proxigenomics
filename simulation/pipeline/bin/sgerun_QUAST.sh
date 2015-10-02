#!/bin/bash

#
# Align query to subject
#

#$ -e logs/
#$ -o logs/
#$ -cwd
#$ -N QUASTJOB

if [ -z "$JOB_ID" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

METAQUAST=external/quast/metaquast.py

if [ -z "$JOB_ID" ] # SUBMIT MODE
then
	if [ $# -ne 3 ]
	then
		echo "Usage: <REFSEQ> <CONTIGS> <OUTPUT REPORT PATH>"
		exit 1
	fi

	echo "Submitting run"
	#trap 'rollback_rm_file $4; exit $?' INT TERM EXIT
	qsub -sync yes -V -pe smp 4 -v REF=$1,CTG=$2,OUTPUT=$3 $0
	#trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"
	echo "Ref $REF"
	echo "Ctg $CTG"

	ODIR=`dirname $OUTPUT`
	$METAQUAST -t $NSLOTS -R $REF -o $ODIR/work $CTG & \
		cp $ORIR/work/combined_quast_output/report.tsv $ODIR & \
		tar czf quast.tar.gz $ODIR/work

fi
