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

export PATH=/work/bin:$PATH
export ITHREE_GIT=$HOME/git/ithree

#
#
# PBS configuration choices can go here.
#
#$ -V
#$ -e logs/
#$ -o logs/
#$ -N A5JOB
#$ -cwd

# Path to the a5 executable
# This might need to be set!
A5EXE=$PROXIHOME/simulation/pipeline/external/a5/bin/a5_pipeline.pl

#
# Check if invocation was via PBS queue.
# 
# If this is not being run from the queue then
# act as a regular shell script and provide a 
# user interface. Successful invocation will
# then cause the script to submit itself to
# the queue.
#
if [ -z "$JOB_ID" ]
then
	####################
	# Regular CLI mode #
	####################

	# defaults
	TAG=a5
	NCPU=4
	MEM=32gb
	QUEUE=all.q
	OVERWRITE="false"

	while getopts ":mt:wfn:M:q:" opt
	do
		case $opt in
			m)
				METAGENOME="--metagenome"
				;;
			w)
				WAITOPT="-sync yes"
				;;
			f)
				OVERWRITE="true"
				;;
			t)
		                TAG=$OPTARG
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
	
	if [ $# -ne 2 ] && [ $# -ne 3 ]
	then
		echo ""
		echo "Submit an A5 job to the queue at UTS in two ways"
		echo ""
		echo "Single set of paired Reads"
		echo "" 
		echo " Usage: [OPTIONS] <READS1> <READS2> <OUTPUT_PATH>"
		echo ""
		echo "Using a libfile"
		echo ""
		echo " Usage: [OPTIONS] <LIBFILE> <OUTPUT_PATH>"
		echo ""
		echo "Optional Arguments"
		echo ""
		echo " -m    enable metagenome behaviour in A5 (optional)"
		echo " -w    block on submission, wait for assembly to finish (optional)"
		echo " -f    force overwriting of output folder contents (optional)"
		echo " -t    prefix tag for output file names"
                echo " -n    number of concurrent threads [default: 2]"
		echo " -M    job memory with units. [default: 32gb]"
		echo " -q    job queue [default: smallq]"
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

		FINAL_FILE_CHECK=`find $3 -maxdepth 1 -type f -name ${TAG}.contigs.fasta ! -empty`
		if [ -n "$FINAL_FILE_CHECK" ]
		then
			echo "Final assembly contigs found in $3, skipping execution"
			exit 0
		else
			echo "Final assembly contigs NOT FOUND in $3, going forward to submission"
		fi
	
		# Queue submission
		CMD=`readlink -f $0`
		qsub -l h=fat-centos -q $QUEUE -pe smp $NCPU $WAITOPT -v NCPU=$NCPU,TAG=$TAG,OVERWRITE=$OVERWRITE,OPTIONS=$METAGENOME,READ1=$R1,READ2=$R2,OUTDIR=$3 $CMD
		
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
		CMD=`readlink -f $0`
		qsub -q $QUEUE -pe smp $NCPU $WAITOPT -v NCPU=$NCPU,TAG=$TAG,OVERWRITE=$OVERWRITE,OPTIONS=$METAGENOME,LIBFILE=$LIBFILE,OUTDIR=$2 $CMD
	fi
else
	##############
	# Queue mode #
	##############
	
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
		$A5EXE --threads=$NCPU $OPTIONS $READ1 $READ2 $TAG
	else
		# Libfile invocation
		$A5EXE --threads=$NCPU $OPTIONS $LIBFILE $TAG
	fi
fi
