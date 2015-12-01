#!/bin/bash

#$ -e logs/
#$ -o logs/
#$ -cwd
#$ -N SPADESJOB

if [ -z "$JOB_ID" ]
then
    BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
    source $BINDIR/bash_init.sh
fi

SPADES_EXE=external/spades/bin/spades.py

if [ -z "$JOB_ID" ] # SUBMIT MODE
then

	# defaults
	NCPU=4
	MEM=32
	QUEUE=all.q
	OVERWRITE="false"

	while getopts ":cwfn:M:q:" opt
	do
		case $opt in
			c)
				CAREFUL="--careful"
				;;
			w)
				WAITOPT="-sync yes"
				;;
			f)
				OVERWRITE="true"
				;;
            n)
                NCPU=$OPTARG
                ;;
			M)
				MEM=$OPTARG
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

    if [ $# -ne 2 ]
    then
		echo "Usage: [options] <R1> <R2> <OUTDIR>"
		echo ""
		echo "Optional Arguments"
		echo ""
		echo " -c    enable careful option in assembly"
		echo " -w    block on submission, wait for assembly to finish (optional)"
		echo " -f    force overwriting of output folder contents (optional)"
        echo " -n    number of concurrent threads [default: 2]"
		echo " -M    job memory in GB without units. [default: 32 GB]"
		echo " -q    job queue [default: $QUEUE]"
		echo ""
		echo " READS1     first input reads file (paired reads)"
		echo " READS2     second input reads file (paired reads)"
		echo " LIBFILEA5  library file for multiple read-sets"
		echo " OUTPUT_PATH containing output folder"
		echo ""
		exit 1
	fi
	echo "Submitting run"

	R1=`readlink -f $1`
	if [ ! -f $R1 ] || [ -z "$R1" ]
	then
		echo "$1 not found"
		exit 1
	fi

	R2=`readlink -f $2`
	if [ ! -f $R2 ] || [ -z "$R2" ]
	then
		echo "$2 not found"
		exit 1
	fi

	FINAL_FILE_CHECK=`find $3 -maxdepth 1 -type f -name contigs.fasta ! -empty`
	if [ -n "$FINAL_FILE_CHECK" ]
	then
		echo "Final assembly contigs found in $3, skipping execution"
		exit 0
	else
		echo "Final assembly contigs NOT FOUND in $3, going forward to submission"
	fi

	# Queue submission
	CMD=`readlink -f $0`
	qsub -q $QUEUE -pe smp $NCPU $WAITOPT -l mem=$MEM -v MEM=$MEM,NCPU=$NCPU,TAG=$TAG,OVERWRITE=$OVERWRITE,OPTIONS=$CAREFUL,READ1=$R1,READ2=$R2,OUTDIR=$3 $CMD
	echo "Finished"

else # EXECUTION MODE
	echo "Running"

	# Preserve pre-existing directories
	if [ "$OVERWRITE" != "true" ] && [ -e $OUTDIR ]
	then
		echo "${OUTDIR}: output containing folder already exists"
		exit 1
	fi
	mkdir -p $OUTDIR

	# Run A5 pipeline
    $SPADES_EXE --threads $NCPU --memory $MEM $OPTIONS -1 $READ1 -2 $READ2 -o $OUTDIR
fi
