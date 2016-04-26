#!/bin/bash

#
# Map sequences in simulation
#

#$ -e logs/
#$ -o logs/
#$ -cwd
#$ -N SRMCLJOB

if [ -z "$JOB_ID" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

CMD=external/srmcl

if [ -z "$JOB_ID" ] # SUBMIT MODE
then
	while getopts ":b:i:w:p:q:t:" opt
	do
		case $opt in
			b)
				CMDOPT="$CMDOPT -b $OPTARG"
				;;
			i)
				CMDOPT="$CMDOPT -i $OPTARG"
				;;
			w)
				CMDOPT="$CMDOPT -w $OPTARG"
				;;
			p)
				CMDOPT="$CMDOPT -p $OPTARG"
				;;
			q)
				CMDOPT="$CMDOPT -q $OPTARG"
				;;
			t)
				CMDOPT="$CMDOPT -t $OPTARG"
				;;
		esac
	done

	shift $(( OPTIND-1 ))

	if [ $# -ne 2 ]
	then
		echo "Usage: [-b balance] [-i inflation] [-w quality] [-p redundancy] [-q penalty] [-t iterations] <input> <output>"
		exit 1
	fi
	
	export CMDOPT
	
	echo "Submitting run"
	#trap 'rollback_rm_file $2; exit $?' INT TERM EXIT
	qsub -sync yes -V -v CMDOPT="$CMDOPT",INPUT=$1,OUTPUT=$2 $0
	#trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"

	if [ -s $INPUT ]
	then
		$CMD $CMDOPT -o $OUTPUT $INPUT
	else
		echo "Input had no data" > $OUTPUT
	fi

fi
