#!/usr/bin/env python
from Bio import SeqIO

import argparse
import os
import subprocess
import sys
import numpy
from numpy import random

TMP_INPUT = 'seq.'+str(numpy.random.randint(0,999999999))+'.tmp'
TMP_OUTPUT = 'reads.'+str(numpy.random.randint(0,999999999))+'.tmp'
R1_FILE = '{0}1.fq'.format(TMP_OUTPUT)
R2_FILE = '{0}2.fq'.format(TMP_OUTPUT)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Simulate a metagenomic data set from an abundance profile')
    parser.add_argument('-M', '--max-coverage', metavar='INT', type=int, required=True,
                        help='Coverage of must abundant taxon')
    parser.add_argument('-S', '--seed', metavar='INT', type=int, required=True, help='Random seed')
    parser.add_argument('-l', '--read-len', metavar='INT', type=int, required=True, help='Read length')
    parser.add_argument('-m', '--insert-len', metavar='INT', type=int, required=True, help='Insert length')
    parser.add_argument('-s', '--insert-sd', metavar='INT', type=int, required=True, help='Insert standard deviation')
    parser.add_argument('--art-path', default='ART_illumina', help='Path to ART executable [default: ART_illumina]')
    parser.add_argument('--log', default='metaART.log', type=argparse.FileType('w'), help='Log file name')
    parser.add_argument('-n', '--num-samples', metavar='INT', type=int, required=True, help='Number of transect samples')
    parser.add_argument('-U', '--lognorm-ra-mu', metavar='FLOAT', type=float, required=True, help='Lognormal relative abundance mu parameter')
    parser.add_argument('-u', '--lognorm-ra-sigma', metavar='FLOAT', type=float, required=True, help='Lognormal relative abundance sigma parameter')
    parser.add_argument('fasta', metavar='MULTIFASTA',
                        help='Input multi-fasta of all sequences')
    parser.add_argument('output_base', metavar='OUTPUT BASE',
                        help='Output file name')
    args = parser.parse_args()


    seq_index = SeqIO.index(args.fasta, 'fasta')
    numpy.random.seed(args.seed)
    all_R1 = open('{0}1.fq'.format(args.output_base), 'w')
    all_R2 = open('{0}2.fq'.format(args.output_base), 'w')

    # generate N simulated communities
    for n in range(0,args.num_samples):

        # generate a lognormal relative abundance profile
        profile = {}
        ra_sum = 0
        for seq_id in seq_index:
            profile[seq_id] = numpy.random.lognormal(args.lognorm_ra_mu, args.lognorm_ra_sigma)
            ra_sum += profile[seq_id]

        for seq_id in seq_index:
            profile[seq_id] /= ra_sum
        print "Sample " + str(n) + " rel abundances " + ", ".join(map(str, profile))

        output_R1_name = '{0}.{1}.r1.fq'.format(args.output_base,n)
        output_R2_name = '{0}.{1}.r2.fq'.format(args.output_base,n)
        with open(output_R1_name, 'w') as output_R1, \
                open(output_R2_name, 'w') as output_R2:
            try:
                for seq_id in profile:
                    coverage = float(profile[seq_id] * args.max_coverage)
                    print 'Requesting {0} coverage for {1}'.format(coverage, seq_id)
                    ref_seq = seq_index[seq_id]
                    ref_len = len(ref_seq)
                    SeqIO.write([ref_seq], TMP_INPUT, 'fasta')
                    subprocess.call([args.art_path,
                                     '-p',   # paired-end sequencing
                                     '-na',  # no alignment file
                                     '-rs', str(args.seed),
                                     '-m', str(args.insert_len),
                                     '-s', str(args.insert_sd),
                                     '-l', str(args.read_len),
                                     '-f', str(coverage),
                                     '-i', TMP_INPUT,
                                     '-o', TMP_OUTPUT], stdout=args.log, stderr=args.log)

                    # count generated reads
                    r1_n = 0
                    for seq in SeqIO.parse(R1_FILE, 'fastq'):
                        r1_n += 1

                    r2_n = 0
                    for seq in SeqIO.parse(R1_FILE, 'fastq'):
                        r2_n += 1

                    effective_cov = args.read_len * (r1_n + r2_n) / float(ref_len)
                    print 'Generated {0} paired-end reads for {1}, {2:.3f} coverage'.format(r1_n, seq_id, effective_cov)
                    if r1_n != r2_n:
                        print 'Error: paired-end counts do not match {0} vs {1}'.format(r1_n, r2_n)
                        sys.exit(1)

                    with open(R1_FILE, 'r') as tmp_h:
                        all_R1.write(tmp_h.read())

                    with open(R2_FILE, 'r') as tmp_h:
                        all_R2.write(tmp_h.read())

                    with open(R1_FILE, 'r') as tmp_h:
                        output_R1.write(tmp_h.read())
                        os.remove(tmp_h.name)

                    with open(R2_FILE, 'r') as tmp_h:
                        output_R2.write(tmp_h.read())
                        os.remove(tmp_h.name)
                    os.remove(TMP_INPUT)

            finally:
                pass
            output_R1.close()
            output_R2.close()
            os.system("gzip " + output_R1_name)
            os.system("gzip " + output_R2_name)

    all_R1.close()
    all_R2.close()
    os.system("gzip " + all_R1.name)
    os.system("gzip " + all_R2.name)
