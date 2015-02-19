#!/bin/bash

if [ -z "$PBS_ENVIRONMENT" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

#
# Align query to subject
#

#PBS -q smallq
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N LASTJOB

LASTAL=$HOME/bin/lastal
LASTDB=$HOME/bin/lastdb
MAFCONV=$HOME/bin/maf-convert.py

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then
	if [ $# -ne 3 ]
	then
		echo "Usage: <subject> <query> <output_base>"
		exit 1
	fi

	echo "Submitting run"
	trap 'rollback_rm_file $3; exit $?' INT TERM EXIT
	qsub -W block=true -v SUBJECT=$1,QUERY=$2,OUTPUT=$3 $0
	trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR

    if [ ! -e ${SUBJECT}.prj ]
    then
        echo "There doesn't appaear to be a LAST DB associated with this reference data"
        echo "Creating database..."
        $LASTDB $FASTA $FASTA
    fi

	# align
	$LASTAL $SUBJECT $QUERY | $MAFCONV psl > $OUTPUT
fi
