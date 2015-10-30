#!/usr/bin/env python
from __future__ import division
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
num_strains = int(sys.argv[1])
ref_fa = sys.argv[2]
num_samples = len(sys.argv) - 4

depths = dict()
for a in alphabet:
    depths[a] = dict()
variant_sites = list()

##
# make variant calls for each input bam
# parse the VCF and store calls in snvs dict
#
for i in range(num_samples):
    cur_vcf = str(i) + ".vcf"
    lofreq_cmd = lofreq + " call " + "--no-default-filter -f " + ref_fa + " -o " + cur_vcf + " " + sys.argv[i+3] 
    print lofreq_cmd
    os.system(lofreq_cmd)
    vcf_file = open(cur_vcf)
    for line in vcf_file:
        if line.startswith("#"):
            continue
        line.rstrip()
        d = line.split("\t")
        if not d[6].startswith("PASS"):
            continue    # didnt pass filters
        locus = d[0] + "\t" + d[1]  # chrom & site
        m = re.search('DP=(.+);AF=(.+);SB', d[7])
        vac = int(float(m.group(1)) * float(m.group(2)))
        if not locus in depths[alphabet[0]]:
            for a in alphabet:
                depths[a][locus] = dict()
                for j in range(num_samples):
                    depths[a][locus][j] = 0
        depths[d[4]][locus][i] = vac
        depths[d[3]][locus][i] = int(m.group(1)) - vac
        variant_sites = variant_sites + [locus]

##
# write out a file with SNVs and sample count for Bayesian PNMF
#
num_sites = len(variant_sites)
snv_filename = "snv_file.data.R"
snv_file = open(snv_filename, "w")
snv_file.write("U<-" + str(num_sites) + "\n")  # number of sites
snv_file.write("T<-" + str(num_samples) + "\n")  # number of time points
snv_file.write("S<-" + str(num_strains) + "\n")  # number of time points
nota = "nota <- c("
notc = "notc <- c("
notg = "notg <- c("
nott = "nott <- c("
sepchar = ""
for site in variant_sites:
    for i in range(num_samples):
        nota = nota + sepchar + str(depths['A'][site][i])
        notc = notc + sepchar + str(depths['C'][site][i])
        notg = notg + sepchar + str(depths['G'][site][i])
        nott = nott + sepchar + str(depths['T'][site][i])
        sepchar = ","

snv_file.write(nota+")\n")
snv_file.write(notc+")\n")
snv_file.write(notg+")\n")
snv_file.write(nott+")\n")
snv_file.close()

##
# run the Poisson NMF
#
bpnmf_filename = "decon.csv"
bpnmf_cmd = "genotypes_acgt variational data file=" + snv_filename + " output file=" + bpnmf_filename
os.system(bpnmf_cmd)

##
# summarize the tip partials and create a BEAST XML
#
beast_filename = "beast.xml"
bpnmf_file = open(bpnmf_filename)
beast_file = open(beast_filename, "w")

# init 2D array of tip partials
tip_partials = [[[0 for x in range(num_sites)] for s in range(num_strains)] for x in range(len(alphabet))]

ll = -2 # skip the first line (csv header)
for line in bpnmf_file:
    if line.startswith("#"):
        continue 
    ll += 1
    if ll < 0:
        continue
    d = line.split(",")
    for i in range(len(alphabet)):
        for s in range(num_strains):
            begin = num_samples + 1 + num_sites * num_strains * i + num_sites * s
            end = begin + num_sites
            for j in range(begin,end):
                tip_partials[i][s][j-begin] += float(d[j])

    # normalize to a tip partial distribution
    for s in range(num_strains):
        for j in range(num_sites):
            m = 0
            for i in range(len(alphabet)):
                m = tip_partials[i][s][j] if tip_partials[i][s][j] > m else m 
            for i in range(len(alphabet)):
                tip_partials[i][s][j] = tip_partials[i][s][j] / m

# write to xml
for s in range(num_strains):
    for i in range(len(alphabet)):
        beast_file.write("<TipPartials id=\"" + str(s) + "\" char=\""+ str(alphabet[i]) + "\">\n")
        beast_file.write(",".join(map(str,tip_partials[i][s])))
        beast_file.write("\n</TipPartials>\n")
   
beast_file.close()

##
# run BEAST
#

