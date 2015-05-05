#!/bin/bash

#
# Create graph files
#

#PBS -q smallq
#PBS -l select=1:ncpus=1:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N GRAPHJOB

if [ -z "$PBS_ENVIRONMENT" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

BAMTOEDGES=bin/bamToEdges.py

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then

    while getopts ":s" opt
	do
		case $opt in
			s)
				SELFLOOP=1
				;;
			\?)
				echo "Invalid option: -$OPTARG"
				;;
		esac
	done

	shift $(( OPTIND-1 ))

	if [ $# -ne 3 ]
	then
		echo "Usage: [options] [hic2ctg.bam] [edge out] [node out]"
		echo "options:"
		echo "  -s   add self-loops for all nodes"
		exit 1
	fi

	echo "Submitting run"
	TARGET=( $2 $3 )
	trap 'rollback_rm_files ${TARGET[@]}; exit $?' INT TERM EXIT
	qsub -W block=true -v SELFLOOP=$SELFLOOP,HIC2CTG=$1,EDGES=$2,NODES=$3 $0
	trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR

    if [ -z $SELFLOOP ]
    then
    	$BAMTOEDGES --add-selfloops $HIC2CTG $EDGES $NODES
    else
    	$BAMTOEDGES --add-selfloops $HIC2CTG $EDGES $NODES
    fi
	
fi
