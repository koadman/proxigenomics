#!/bin/bash

#
# Create graph files
#

#$ -e logs/
#$ -o logs/
#$ -cwd
#$ -N GRAPHJOB

if [ -z "$JOB_ID" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

PROXIHOME=~/git/proxigenomics/simulation/pipeline
BAMTOEDGES=$PROXIHOME/bin/bamToEdges.py

if [ -z "$JOB_ID" ] # SUBMIT MODE
then
    MINLEN=1000
    WAITOPT="-sync no"

    while getopts ":swm:" opt
	do
		case $opt in
            w)
                WAITOPT="-sync yes"
                ;;
			s)
				SELFLOOP=1
				;;
            m)
                MINLEN=$OPTARG
                ;;
			\?)
				echo "Invalid option: -$OPTARG"
				;;
		esac
	done

	shift $(( OPTIND-1 ))

	if [ $# -ne 4 ]
	then
		echo "Usage: [options] [hic2ctg.bam] [edge out] [node out] [graph out]"
		echo "options:"
        echo "  -w          wait until job is completed"
		echo "  -s          add self-loops for all nodes"
        echo "  -m LENGTH   minimum contig length to include in graph [1000]"
		exit 1
	fi

	echo "Submitting run"
	qsub $WAITOPT -V -v SELFLOOP=$SELFLOOP,HIC2CTG=$1,MINLEN=$MINLEN,EDGES=$2,NODES=$3,GRAPH=$4 $0
	echo "Finished"

else # EXECUTION MODE
	echo "Running"

    if [ -z $SELFLOOP ]
    then
    	$BAMTOEDGES --add-selfloops --minlen $MINLEN --graphml $GRAPH $HIC2CTG $EDGES $NODES
    else
    	$BAMTOEDGES --add-selfloops --minlen $MINLEN --graphml $GRAPH $HIC2CTG $EDGES $NODES
    fi

fi
