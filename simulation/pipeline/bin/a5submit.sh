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
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -q smallq

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
	
	while getopts ":mtwf" opt
	do
		case $opt in
			m)
				METAGENOME="--metagenome"
				;;
			t)
				TIMESTAMP="true"
				;;
			w)
				WAITOPT="-W block=true"
				;;
			f)
				OVERWRITE="true"
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
		echo " Usage: [-m] <READS1> <READS2> <OUTPUT_TAG>"
		echo ""
		echo "Using a libfile"
		echo ""
		echo " Usage: [-m] <LIBFILE> <OUTPUT_TAG>"
		echo ""
		echo "Arguments"
		echo ""
		echo " -m    enable metagenome behaviour in A5 (optional)"
		echo " -t    append a timestamp to the output folder (optional)"
		echo " -w    block on submission, wait for assembly to finish (optional)"
		echo " -f    force overwriting of output folder contents (optional)"
		echo ""
		echo " READS1     first input reads file (paired reads)"
		echo " READS2     second input reads file (paired reads)"
		echo " LIBFILEA5  library file for multiple read-sets"
		echo " OUTPUT_TAG prefix tag used in output file names (this is not a path)"
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
	
		TESTBASE=`basename $3`
		if [ ${#3} -ne ${#TESTBASE} ]
		then
			echo "$3 should not be a path, please remove any slashes (/)"
			exit 1
		fi

		# Queue submission
		qsub $WAITOPT -N A5job -v TIMESTAMP=$TIMESTAMP,OPTIONS=$METAGENOME,READ1=$R1,READ2=$R2,OUTBASE=$3 $0
		
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
		qsub $WAITOPT -N A5job -v TIMESTAMP=$TIMESTAMP,OPTIONS=$METAGENOME,LIBFILE=$LIBFILE,OUTBASE=$2 $0
	fi
else
	##############
	# Queue mode #
	##############
	
	# Change to the working directory where the job was submitted to the queue
	cd $PBS_O_WORKDIR

	# Create an output directory to contain the results
	if [ "$TIMESTAMP" == "true" ]
	then
		OUTDIR=${OUTBASE}"_"`date +%Y%m%d.%S`
	else
		
		OUTDIR=${OUTBASE}
	fi
	
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
		$A5EXE $OPTIONS $READ1 $READ2 $OUTBASE
	else
		# Libfile invocation
		$A5EXE $OPTIONS $LIBFILE $OUTBASE
	fi
fi
