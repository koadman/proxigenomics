#!/bin/bash
#
# Simple script to submit A5 jobs to PBS queue at UTS
#
# MZD 2014/02/03
#
# This script expects 3 variables passed in by the -v command.
# READ1 - R1 fastq
# READ2 - R2 fastq
# OUTBASE - Base directory for output, also passed to A5 as "base" command
#
# Work will be performed relative to the submission current directory. Paths 
# to data reads and output can be expressed locally.

#
#
# PBS configuration choices can go here.
#
#PBS -q smallq
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N A5JOB

# Path to the a5 executable
# This might need to be set!
A5EXE=$ITHREE_GIT/software/a5_miseq_linux/bin/a5_pipeline.pl


#
# Check if invocation was via PBS queue.
# 
# If this is not being run from the queue then
# act as a regular shell script and provide a 
# user interface. Successful invocation will
# then cause the script to submit itself to
# the queue.
#
if [ -z "$PBS_ENVIRONMENT" ]
then
	####################
	# Regular CLI mode #
	####################

    # defaults
    TAG=a5

	while getopts ":mt:wf" opt
	do
		case $opt in
			m)
				METAGENOME="--metagenome"
				;;
			w)
				WAITOPT="-W block=true"
				;;
			f)
				OVERWRITE="true"
				;;
			t)
		        TAG=$OPTARG
		        ;;
			\?)
				echo "Invalid option: -$OPTARG"
				;;
		esac
	done
	
	shift $(( OPTIND-1 ))
	
	if [ $# -ne 2 ] && [ $# -ne 3 ]
	then
		echo ""
		echo "Submit an A5 job to the queue at UTS in two ways"
		echo ""
		echo "Single set of paired Reads"
		echo "" 
		echo " Usage: [-m] <READS1> <READS2> <OUTPUT_PATH>"
		echo ""
		echo "Using a libfile"
		echo ""
		echo " Usage: [-m] <LIBFILE> <OUTPUT_PATH>"
		echo ""
		echo "Arguments"
		echo ""
		echo " -m    enable metagenome behaviour in A5 (optional)"
		echo " -w    block on submission, wait for assembly to finish (optional)"
		echo " -f    force overwriting of output folder contents (optional)"
		echo " -t    prefix tag for output file names"
		echo ""
		echo " READS1     first input reads file (paired reads)"
		echo " READS2     second input reads file (paired reads)"
		echo " LIBFILEA5  library file for multiple read-sets"
		echo " OUTPUT_PATH containing output folder"
		echo ""
		exit 1
	fi
	
	# Set up paths for read-set submission
	if [ $# -eq 3 ]
	then
		# Check for input read existence
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
	
		# Queue submission
		qsub $WAITOPT -v TAG=$TAG,OVERWRITE=$OVERWRITE,OPTIONS=$METAGENOME,READ1=$R1,READ2=$R2,OUTDIR=$3 $0
		
	# Set up paths for libfile submission
	else
		LIBFILE=`readlink -f $1`
		if [ ! -f $LIBFILE ]
		then
			echo "$1: libfile not found"
			exit 1
		fi
	
		TESTBASE=`basename $2`
		if [ ${#2} -ne ${#TESTBASE} ]
		then
			echo "$2 should not be a path, please remove any slashes (/)"
			exit 1
		fi

		# Queue submission
		qsub $WAITOPT -v TAG=$TAG,OVERWRITE=$OVERWRITE,OPTIONS=$METAGENOME,LIBFILE=$LIBFILE,OUTDIR=$2 $0
	fi
else
	##############
	# Queue mode #
	##############
	
	# Change to the working directory where the job was submitted to the queue
	cd $PBS_O_WORKDIR

	# Preserve pre-existing directories
	if [ "$OVERWRITE" != "true" ] && [ -e $OUTDIR ]
	then
		echo "${OUTDIR}: output containing folder already exists"
		exit 1
	fi
	mkdir -p $OUTDIR
	cd $OUTDIR

	# Run A5 pipeline
	if [ -z "$LIBFILE" ]
	then
		# Read-set invocation
		$A5EXE --threads=2 $OPTIONS $READ1 $READ2 $TAG
	else
		# Libfile invocation
		$A5EXE --threads=2 $OPTIONS $LIBFILE $TAG
	fi
fi
