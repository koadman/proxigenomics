#!/bin/bash

if [ -z "$PBS_ENVIRONMENT" ]
then
	BINDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	source $BINDIR/bash_init.sh
fi

#
# Map sequences in simulation
#


#PBS -q smallq
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -e logs/
#PBS -o logs/
#PBS -N BWAJOB

BWAEXE=$HOME/bin/bwa-0.7.6a/bwa

if [ -z "$PBS_ENVIRONMENT" ] # SUBMIT MODE
then
	if [ $# -ne 6 ]
	then
		echo "Usage: <ref seq> <wgs base> <hic base> <scf2ref> <wgs2scf> <hic2scf>"
		exit 1
	fi
	echo "Submitting run"
	qsub -W block=true -v REFSEQ=$1,WGS_BASE=$2,HIC_BASE=$3,SCF2REF=$4,WGS2SCF=$5,HIC2SCF=$6 $0

else # EXECUTION MODE
	echo "Running"
	cd $PBS_O_WORKDIR
	
	# cp reference sequences to working dir
	cp $REFSEQ .
	LOCALREF=`basename $REFSEQ`

	# just find the final scaffold file from A5
	SCAFFOLDS=`find ${WGS_BASE} -type f -name 'wgs.contigs.fasta'`
	#SCAFFOLDS=`find ${WGS_BASE} -type f -name 'wgs.final.scaffolds.fasta'`

	# create indexes
	$BWAEXE index $LOCALREF
	$BWAEXE index $SCAFFOLDS

	# map
	$BWAEXE mem -t 2 $LOCALREF $SCAFFOLDS > $SCF2REF.sam
	$BWAEXE mem -t 2 $SCAFFOLDS ${WGS_BASE}1.fq ${WGS_BASE}2.fq > $WGS2SCF.sam
	$BWAEXE mem -t 2 $SCAFFOLDS ${HIC_BASE}.fasta > $HIC2SCF.sam
fi
