#!/bin/bash

#
# Align query to subject
#

#PBS -q smallq
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N LASTJOB

if [ -z "$PBS_ENVIRONMENT" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

LASTAL=$ITHREE_GIT/software/bin/lastal
LASTDB=$ITHREE_GIT/software/bin/lastdb
LASTMAP=$ITHREE_GIT/software/bin/last-map-probs
MAFCONV=$ITHREE_GIT/software/bin/maf-convert.py

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then
	if [ $# -ne 3 ]
	then
		echo "Usage: <subject> <query> <output_base>"
		exit 1
	fi

	echo "Submitting run"
	TARGET=( ${1}.prj ${3} )
	trap 'rollback_rm_files ${TARGET[@]}; exit $?' INT TERM EXIT
	qsub -W block=true -v SUBJECT=$1,QUERY=$2,OUTPUT=$3 $0
	trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR

    if [ ! -e ${SUBJECT}.mapdb.prj ]
    then
        echo "There doesn't appaear to be a LAST DB associated with this reference data"
        echo "Creating database..."
        $LASTDB -m1111110 ${SUBJECT}.mapdb $SUBJECT
    fi

	# align
	$LASTAL -Q1 -e120 ${SUBJECT}.mapdb $QUERY | $LASTMAP | $MAFCONV psl > $OUTPUT
fi
