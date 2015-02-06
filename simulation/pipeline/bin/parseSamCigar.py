#!/usr/bin/env python
from Bio import SeqIO
import re
import sys

matcher = re.compile(r'([0-9]+)M')


def count_aligned(cigar):
    return sum([int(s) for s in matcher.findall(cigar)])


class Alignment:
    def __init__(self, bases=0):
        self.bases = bases

    def __repr__(self):
        return repr((self.bases))

    def __str__(self):
        return str(self.bases)

    def add_bases(self, bases):
        self.bases += bases


if len(sys.argv) != 5:
    print 'Usage: [min length] [query fasta] [sam file] [outfile]'
    sys.exit(1)

# Set minimum sequence length and calculate the length of all input fasta sequence
min_length = int(sys.argv[1])
seq_length = {rec.id: len(rec) for rec in SeqIO.parse(sys.argv[2], 'fasta')}

# Parse the SAM file for alignment lengths
taxons = {}
i = 0
with open(sys.argv[3], 'r') as h_in:
    for l in h_in:

        if l.startswith('@'):
            continue

        fields = l.rsplit()
        tx = taxons.get(fields[2])
        if tx is None:
            tx = {}
            taxons[fields[2]] = tx
        txal = tx.get(fields[0])
        if txal is None:
            txal = Alignment()
            tx[fields[0]] = txal
        txal.add_bases(count_aligned(fields[5]))

#
# We need to decide on assignment.
#
# There are contigs which are assigned to more than one taxon.
# What shall be done with these? Simply taking the largest alignmetn
# as the winner might be acceptable for large differences, but I've
# seen examples where both are significant alignments and picking one
# is quite arbitrary.
#

# PICK THE WINNER
# For each scaffold, determine the best alignment by length. The alignment subject
# will then be considered the true source.
best = {}
for t, s in taxons.iteritems():
    for sn, aln in s.iteritems():
        bi = best.get(sn)
        if bi is None or aln.bases > bi['aln'].bases:
            best[sn] = {'tx': t, 'aln': aln, 'slen': seq_length[sn]}

with open(sys.argv[4], 'w') as h_out:
    # Write out the list of winners, with additional columns for validation
    for k, v in best.iteritems():
        if v['aln'].bases > min_length:
            h_out.write('{seq_name} {tax_name} {cov:.4}\n'.format(
                seq_name=k, tax_name=v['tx'], cov=float(v['aln'].bases) / float(v['slen'])))
