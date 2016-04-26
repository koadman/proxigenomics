#!/bin/bash

#
# Map sequences in simulation
#

#$ -e logs/
#$ -o logs/
#$ -cwd
#$ -N MAPJOB

if [ -z "$JOB_ID" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

BWA=external/bwa
SAMTOOLS=external/samtools

if [ -z "$JOB_ID" ] # SUBMIT MODE
then

    WAITOPT="-sync no"
	NCPU=4
	QUEUE=all.q

	while getopts ":wn:q:" opt
	do
		case $opt in
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

	if [ $# -ne 3 ] && [ $# -ne 4 ]
	then
		echo "Usage: [options] <subject> <query> <output>"
		echo "or"
		echo "Usage: [options] <subject> <query1> <query2> <output>"
        echo ""
        echo "Options:"
        echo "-w    wait for job to finish"
        echo "-n    number of cpus"
        echo "-q    queue name [all.q]"
		exit 1
	fi

	# get full path of script
	CMD=`readlink -f $0`

	if [ $# -eq 3 ]
	then
    	echo "Submitting single query run"
    	#BASE=${3%.bam}
    	#TARGET=( ${BASE}.sam ${BASE}.bam ${BASE}.bai ${BASE}.idxstats ${BASE}.flagstats )
    	#trap 'rollback_rm_files ${TARGET[@]}; exit $?' INT TERM EXIT
	    qsub $WAITOPT -V -pe smp $NCPU -v SUBJECT=$1,QUERY1=$2,QUERY2="",OUTPUT=$3 $CMD
    	trap - INT TERM EXIT
    elif [ $# -eq 4 ]
    then
    	echo "Submitting double query run"
    	#BASE=${4%.bam}
    	#TARGET=( ${BASE}.sam ${BASE}.bam ${BASE}.bai ${BASE}.idxstats ${BASE}.flagstats )
        #trap 'rollback_rm_files ${TARGET[@]}; exit $?' INT TERM EXIT
        qsub  $WAITOPT -V -pe smp $NCPU -v SUBJECT=$1,QUERY1=$2,QUERY2=$3,OUTPUT=$4 $CMD
    	trap - INT TERM EXIT
    fi
	echo "Finished"

else # EXECUTION MODE
	echo "Running"

    BASE=${OUTPUT%.bam}

	# Map reads. The second query file can be empty
#	$BWA mem -t $NSLOTS $SUBJECT $QUERY1 $QUERY2 > ${BASE}.sam

#    $SAMTOOLS view -@ $NSLOTS -bS ${BASE}.sam | $SAMTOOLS sort -@ $NSLOTS - $BASE

	$BWA mem -t $NSLOTS $SUBJECT $QUERY1 $QUERY2 | \
        $SAMTOOLS view -@ $NSLOTS -bS - | $SAMTOOLS sort -@ $NSLOTS - $BASE

    $SAMTOOLS index $OUTPUT
    $SAMTOOLS idxstats $OUTPUT > ${BASE}.idxstats
    $SAMTOOLS flagstat $OUTPUT > ${BASE}.flagstat

fi
