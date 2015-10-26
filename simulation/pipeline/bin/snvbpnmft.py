#!/usr/bin/env python
#
# wrapper script to: 
# (1) compute SNVs using lofreq*
# (2) compute a variational Bayes Poisson non-negative matrix factorization
# (3) do phylogenetic inference on the strain genotype matrix
#

import sys
import re
import os

lofreq = "lofreq"
alphabet = ['A','C','G','T']

# parse the command-line
if len(sys.argv)<5:
    print "Usage: snvbpnmft.py <number of strains> <reference fasta> <sample 1 bam> <sample 2 bam> .. [sample N bam]";
    sys.exit(-1)
num_strains = sys.argv[0]
ref_fa = sys.argv[1]
num_samples = len(sys.argv) - 2

depths = dict()
for a in alphabet:
    depths[a] = dict()
variant_sites = list()

##
# make variant calls for each input bam
# parse the VCF and store calls in snvs dict
#
for i in range(1,num_samples):
    cur_vcf = str(i) + ".vcf"
    lofreq_cmd = lofreq + " call " + "-r " + ref_fa + " -o " + cur_vcf + " " + sys.argv[i] 
    os.system(lofreq_cmd)
    vcf_file = open(cur_vcf)
    for line in vcf_file:
        if line.startsWith("#"):
            continue
        line.rstrip()
        d = line.split("\t")
        if not d[6].startsWith("PASS"):
            continue    # didnt pass filters
        locus = d[0] + "\t" + d[1]  # chrom & site
        m = re.search('DF=(.+)\;AF=(.+);SB', d[7])
        vac = m.group(0) * m.group(1)
        depths[d[4]][locus][i] = vac
        depths[d[3]][locus][i] = m.group(0) - vac
        variant_sites.add(locus)        

##
# write out a file with SNVs and sample count for Bayesian PNMF
#
snv_filename = "snv_file.data.R"
snv_file = open(snv_filename, "w")
snv_file.write("U<-" + str(len(variant_sites)))  # number of sites
snv_file.write("T<-" + str(num_samples))  # number of time points
snv_file.write("S<-" + str(num_strains))  # number of time points
nota = "nota <- c("
notc = "notc <- c("
notg = "notg <- c("
nott = "nott <- c("
sepchar = ""
for site in variant_sites:
    for i in range(1,num_samples):
        nota = nota + sepchar
        notc = notc + sepchar
        notg = notg + sepchar
        nott = nott + sepchar
        if exists(snv_depths[site][i]):
            nota = nota + str(depths['A'][site][i])
            notc = notc + str(depths['C'][site][i])
            notg = notg + str(depths['G'][site][i])
            nott = nott + str(depths['T'][site][i])
        sepchar = ","

snv_file.write(nota)
snv_file.write(notc)
snv_file.write(notg)
snv_file.write(nott)
snv_file.close()

##
# run the Poisson NMF
#
bpnmf_filename = "decon.csv"
bpnmf_cmd = "genotypes variational data file=" + snv_filename + " output file=" + bpnmf_filename
os.system(bpnmf_cmd)

##
# summarize the tip partials and create a BEAST XML
#
bpnmf_file = open(bpnmf_filename)
tip_partials = dict()
for line in bpnmf_file:
    d = line.split(",")
    
    for i in len(alphabet):
        begin = num_samples + 1 + num_sites * i
        end = begin + num_sites
        for j in range(begin,end):
            tip_partials[ alphabet[i] ] += float(d[j])


