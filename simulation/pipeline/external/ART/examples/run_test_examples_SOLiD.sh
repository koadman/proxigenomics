#!/bin/bash
#SOLiD test examples

art_SOLiD=../art_SOLiD
# 1) simulation of SOLiD 32bp single-end reads at 10X coverage using the built-in SOLiD FLX read  profile

$art_SOLiD -s ./testSeq.fa ./single_SOLiD_test 32 10

# 2) simulation of SOLiD 25bp paried-end reads at 5X coverage with the mean fragment size 500 
#    and standard deviation 20 using the built-in SOLiD FLX read  profile

$art_SOLiD -s ./testSeq.fa ./paired_SOLiD_test 25 5 500 20



