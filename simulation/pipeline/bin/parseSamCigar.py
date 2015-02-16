#!/usr/bin/env python
from Bio import SeqIO
import re
import argparse

# CIGAR matcher for aligning regions
matcher = re.compile(r'([0-9]+)M')


def count_aligned(cigar):
    """
    Sum the length of aligned regions listed in CIGAR pattern
    :param cigar: SAMtools CIGAR string
    :return: total number of aligned bases
    """
    return sum([int(seg_len) for seg_len in matcher.findall(cigar)])


class Alignment:
    """
    Simple class for holding the cumulative alignment length
    of a query sequence on a reference.
    """
    def __init__(self, bases=0):
        self.bases = bases

    def __repr__(self):
        return repr(self.bases)

    def __str__(self):
        return str(self.bases)

    def add_bases(self, bases):
        self.bases += bases


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Determine the location of query sequences on the reference')
    parser.add_argument('--degen', default=False, action='store_true',
                        help='Allow degenerate query assignment')
    parser.add_argument('--min-length', type=int, required=True, help='Minimum length in bp')
    parser.add_argument('query_fasta', metavar='FASTA', nargs=1, help='Query fasta sequences')
    parser.add_argument('sam_file', metavar='SAMFILE', nargs=1,
                        help='SAM file with query sequences aligned to reference')
    parser.add_argument('output_file', metavar='OUTPUT', nargs=1, help='Output file name')
    args = parser.parse_args()

#    if len(sys.argv) != 5:
#        print 'Usage: [min length] [query fasta] [sam file] [outfile]'
#        sys.exit(1)

    # Calculate the length of query sequences
    seq_length = {rec.id: len(rec) for rec in SeqIO.parse(args.query_fasta[0], 'fasta')}

    # Parse the SAM file for alignment lengths, where each
    # query can accumulate total length there are multiple
    # record against the same reference.

    refs_repository = {}  # repo of references, indexed by name
    i = 0
    with open(args.sam_file[0], 'r') as h_in:
        for l in h_in:
            # skip header fields
            if l.startswith('@'):
                continue
            fields = l.rsplit()

            # fetch the reference record from repo
            # if it doesn't exist, add it to repo.
            ref_name = refs_repository.get(fields[2])
            if ref_name is None:
                ref_name = {}
                refs_repository[fields[2]] = ref_name

            # fetch query to ref alignment by qname
            # from an individual reference repo.
            # if it doesn't exist, add alignment to repo.
            query_aln = ref_name.get(fields[0])
            if query_aln is None:
                query_aln = Alignment()
                ref_name[fields[0]] = query_aln
            query_aln.add_bases(count_aligned(fields[5]))

    #
    # We need to decide on assignment.
    #
    # There are query sequences which are assigned to more than one reference.
    # What shall be done with these? Simply taking the largest alignment
    # as the winner might be acceptable for large differences, but I've
    # seen examples where both are significant alignments and picking one
    # is quite arbitrary.
    #

    if args.degen:
        '''
        Assign query sequence to all significant alignments.
        '''

        # TODO

    else:
        '''
        Pick a winner

        For each scaffold, determine the best alignment by length. The alignment subject
        will then be considered the true source.
        '''

        best = {}
        # for all references
        for ref_name, queries in refs_repository.iteritems():
            # for all queries aligning to each reference
            for q_name, aln in queries.iteritems():
                # if first time we've seen query or
                # its longer than the current leader
                # it becomes the new leader
                bi = best.get(q_name)
                if bi is None or aln.bases > bi['qry'].bases:
                    best[q_name] = {'ref': ref_name, 'qry': aln, 'qlen': seq_length[q_name]}

        with open(args.output_file[0], 'w') as h_out:
            # Write out the list of winners, with additional columns for validation
            for k, v in best.iteritems():
                if v['qry'].bases > args.min_length:
                    h_out.write('{seq_name} {tax_name} {cov:.4}\n'.format(
                        seq_name=k, tax_name=v['ref'], cov=float(v['qry'].bases) / float(v['qlen'])))
