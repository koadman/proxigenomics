#!/bin/bash

#
# Align query to subject
#

#$ -e logs/
#$ -o logs/
#$ -cwd
#$ -N LASTJOB

if [ -z "$JOB_ID" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

LASTAL=external/last/bin/lastal
LASTDB=external/last/bin/lastdb
MAFCONV=external/last/bin/maf-convert

if [ -z "$JOB_ID" ] # SUBMIT MODE
then

    WAITOPT="-sync no"
	NCPU=4
	QUEUE=all.q

	while getopts ":rwn:q:" opt
	do
		case $opt in
            r)
                REBUILD="true"
                ;;
			w)
				WAITOPT="-sync yes"
				;;
            n) 
                NCPU=$OPTARG
                ;;
			q)
				QUEUE=$OPTARG
				;;
			\?)
				echo "Invalid option: -$OPTARG"
				;;
		esac
	done
	
	shift $(( OPTIND-1 ))

	if [ $# -ne 3 ]
	then
		echo "Usage: [options] <subject> <query> <output_base>"
        echo ""
        echo "Options:"
        echo "-r    force rebuild of index"
        echo "-w    wait for job to finish"
        echo "-n    number of cpus"
        echo "-q    queue name [all.q]"
		exit 1
	fi


    if [ "$REBUILD" == "true" ]
    then
        rm -f ${SUBJECT}.prj
    fi

	echo "Submitting run"
	#TARGET=( ${1}.prj ${3} )
	#trap 'rollback_rm_files ${TARGET[@]}; exit $?' INT TERM EXIT
	CMD=`readlink -f $0`
	qsub $WAITOPT -q $QUEUE -pe smp $NCPU -V -v SUBJECT=$1,QUERY=$2,OUTPUT=$3 $CMD
	#trap - INT TERM EXIT
	echo "Finished"

else # EXECUTION MODE
	echo "Running"
    
    if [ ! -e ${SUBJECT}.prj ]
    then
        echo "There doesn't appaear to be a LAST DB associated with this reference data"
        echo "Creating database..."
        $LASTDB $SUBJECT $SUBJECT
    fi

	# align
	$LASTAL -P $NSLOTS $SUBJECT $QUERY | $MAFCONV psl > $OUTPUT
fi
