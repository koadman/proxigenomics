#!/bin/bash

if [ $# -ne 2 ]
then
	echo "Usage: [min length] [run folder]"
	exit 1
fi

if [ ! -d $2 ]
then
	echo "$2: simulation run folder not found"
	exit 1
fi

PYTHON=python
#PYTHON=~/python/bin/python

# Produce a ground truth estimation.
# count aligned bases from SAM, then pick the longest alignment as winner for each scaffold in SAM
$PYTHON bin/parseSamCigar.py $1 $2/wgs/wgs.contigs.fasta $2/scf2ref.sam | sort -k1 > truth.txt

# Prepare MCL input
$PYTHON bin/makeMCLinput.py $1 $2

# Run over both weighted and unweighted edge lists
for fn in mclIn.weighted mclIn.unweighted;
do
	start=2.0
	step=0.025
	# step over range in inflation
	for ((i=1; i<=40; i++))
	do
		# make a string representation of floating point inflation parameter in bash
		ip=`echo "scale=4; $start + $step*$i" | bc -l`
		
		# cluster
		mcl $fn --abc -I $ip -o $fn".out."$1"."$ip
		
		# convert mcl output to cluster assignment list and then join with ground truth estimated above
		# we drop assignments which have no ground truth prediction
		awk '{++nc; for (i=1; i<=NF; i++) print $i,"cl"nc;}' $fn".out."$1"."$ip | \
			awk 'BEGIN{while (getline < "truth.txt") truth[$1]=$2; print "scf truth predict"} {if (truth[$1]!="") print $1,truth[$1],$2;}' > $fn".out."$1"."$ip".classification"
	done
done
