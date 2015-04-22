#!/bin/bash

#PBS -q smallq
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N OCLUSTRJOB

if [ -z "$PBS_ENVIRONMENT" ]
then
    BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
    source $BINDIR/bash_init.sh
fi

OCLUSTR=bin/oclustr.py

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then

    # by default, no isolates
    ISOLATES="--no-isolates"

	while getopts ":i" opt
	do
		case $opt in
			i)
				ISOLATES=""
				;;
			\?)
				echo "Invalid option: -$OPTARG"
				;;
		esac
	done

	shift $(( OPTIND-1 ))

    if [ $# -ne 2 ]
    then
		echo "Usage: [-i] <input> <output>"
		echo ""
		echo "-i     include isolated nodes"
		echo ""
		exit 1
	fi
	echo "Submitting run"

	trap 'rollback_rm_file $3; exit $?' INT TERM EXIT
	qsub -W block=true -v OPTIONS=$ISOLATES,INPUT=$1,OUTPUT=$2 $0
	trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR

	if [ -s $INPUT ]
	then
		$OCLUSTR -f mcl $OPTIONS $INPUT $OUTPUT
	else
		echo "Input had no data" > $OUTPUT
	fi

fi
