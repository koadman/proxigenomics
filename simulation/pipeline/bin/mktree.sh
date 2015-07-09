#!/bin/bash

if [ $# -ne 2 ]
then
	echo "Usage: [input multifasta] [output dir]"
	exit 1
fi

sleep 10

INPUT_FASTA=$1
OUTPUT_DIR=$2 

if [ ! -e $OUTPUT_DIR ]
then
	mkdir -p $OUTPUT_DIR
fi

# split into individual files
mkdir $OUTPUT_DIR/split
csplit -f $OUTPUT_DIR/split/seq -b "%0d.fasta" $INPUT_FASTA '%^>%' '/^>/' '{*}' -s

# progressive mauve
progressiveMauve --output=$OUTPUT_DIR/genomes.xmfa $OUTPUT_DIR/split/seq*

# convert to Fasta
convertXmfa.pl -i $OUTPUT_DIR/genomes.xmfa -o $OUTPUT_DIR/genomes.xmfa.fasta -f fasta -c -l NNNNNNNNNN -g xmfa

# Tree
FastTree_accu -nt -gtr $OUTPUT_DIR/genomes.xmfa.fasta > $OUTPUT_DIR/genomes.xmfa.tre

# ANIb
calculate_ani.py --force -m ANIb --threads=2 -i $OUTPUT_DIR/split -o $OUTPUT_DIR/ani
