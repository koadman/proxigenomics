#!/bin/bash
#illumina test examples

# 1) simulation of single-end reads of 35bp with 10X using the built-in combined quality profile

../art_illumina -i ./testSeq.fa -o ./single_end_com -l 35 -f 10 -sam

#convert an aln file to a bed file
../aln2bed.pl single_end_com.bed single_end_com.aln

# 2) simulation of single-end reads of 35bp with 10X using the built-in seperated quality profiles for A, C, G, and T 

../art_illumina -i ./testSeq.fa -o ./single_end_sep -l 35 -f 10 -sp -sam
#convert an aln file to a bed file
../aln2bed.pl single_end_sep.bed single_end_sep.aln

# 3) simulation of paired-end reads of 50bp with the mean fragment size 500 and standard deviation 10
#    using the built-in combined read quality profiles

../art_illumina -i ./testSeq.fa -o ./paired_end_com -l 50 -f 10 -p -m 500 -s 10 -sam
#convert both aln files to a bed file
../aln2bed.pl paired_end_com.bed paired_end_com1.aln paired_end_com2.aln

# 4) simulation of paired-end reads of 50bp with the mean fragment size 500 and standard deviation 10
#   using the built-in seperated quality profiles for A, C, G, and T 

../art_illumina -i ./testSeq.fa -o ./paired_end_sep -l 50 -f 10 -p -m 500 -s 10 -sp -sam
#convert both aln files to a bed file
../aln2bed.pl paired_end_sep.bed paired_end_sep1.aln paired_end_sep2.aln

# 5) simulation of mate-pair reads of 50bp with the mean fragment size 2500 and standard deviation 50
#    using the built-in combined read quality profiles

../art_illumina -i ./testSeq.fa -o ./matepair_com -l 50 -f 10 -p -m 2500 -s 50 -sam

#convert both aln files to a bed file
../aln2bed.pl matepair_com.bed matepair_com1.aln matepair_com2.aln

