#!/usr/bin/env python
from Bio import SeqIO
from collections import OrderedDict

import re
import argparse
import truthtable as tt


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
    Object representing an alignment taken from some alignment format file.

    This object's uses value identity, governed by the values of query and reference names
    NOT memory reference.

    Not all files have direct access to all parameters. In particular, unless
    specified the alignment will not have percentage identity. This will be
    avoided if not defined.
    """
    def __init__(self, query_name, ref_name, bases, query_length, perid=None):
        self.query_name = query_name
        self.ref_name = ref_name
        if not isinstance(bases, (int, long)):
            raise TypeError('bases parameter must be an integer')
        self.align_length = bases
        self.query_length = query_length
        self.perid = perid

    def __repr__(self):
        return self.query_name + ' ' + self.ref_name

    def __str__(self):
        if self.perid:
            return '{0} {1} {2:.4} {3} {4} {5:.4}'.format(
                self.query_name, self.ref_name, self.coverage,
                self.query_length, self.align_length, self.perid)
        else:
            return '{0} {1} {2:.4} {3} {4}'.format(
                self.query_name, self.ref_name, self.coverage,
                self.query_length, self.align_length)


    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return self.query_name == other.query_name and self.ref_name == other.ref_name

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(repr(self))

    def add_bases(self, bases):
        self.align_length += bases

    @property
    def coverage(self):
        return float(self.align_length) / float(self.query_length)


#def ordereddict_rep(dumper, data):
#    return dumper.represent_mapping(u'tag:yaml.org,2002:map', data.items(), flow_style=False)


def write_output(simple, fmt, minlen, mincov, minid, alignment_list, file_name):
    """
    Write result to file in a particular format.

    :param simple: use simple hard tt
    :param fmt: format type
    :param minlen: minimum query length
    :param mincov: minimum alignment coverage
    :param minid: minimum percentage identity
    :param alignment_list: alignment list
    :param file_name: output file name
    """
    if fmt == 'flat':
        write_alignments(minlen, mincov, minid, alignment_list, file_name)
    elif fmt == 'yaml':
        write_truth(simple, minlen, mincov, alignment_list, file_name)


def write_truth(simple, minlen, mincov, alignment_list, file_name):
    """
    Write YAML format truth table, where all relations between queries and subjects
    are included on a single line.

    Eg.

    seq1: [ctgA, ctgB, ctgC]
    seq2: [ctg111, ctg1]
    seq3: [ctg99]

    :param minlen: minimum query length
    :param mincov: minimum alignment coverage
    :param minid: minimum percentage identity
    :param alignment_list: alignment list
    :param file_name: output file name
    """
    truth = tt.TruthTable()
    for aln in alignment_list:
        if aln.align_length > minlen and aln.coverage > mincov:
            t = truth.get(aln.query_name)
            if t is None:
                t = {}
                truth.put(aln.query_name, t)
            t[aln.ref_name] = aln.align_length

    if simple:
        truth.write_hard(file_name)
    else:
        truth.write(file_name)


def write_alignments(minlen, mincov, minid, alignment_list, file_name):
    """
    Write a simple list of assignments, one relation per line.

    Eg.
    seq1 ctgA
    seq1 ctgB
    seq2 ctg111

    :param minlen: minimum query length
    :param mincov: minimum alignment coverage
    :param minid: minimum percentage identity
    :param alignment_list: alignment list
    :param file_name: output file name
    """
    with open(args.output_file[0], 'w') as h_out:
        for aln in alignment_list:
            if aln.align_length > minlen and aln.coverage > mincov:
                h_out.write(str(aln) + '\n')
            else:
                h_out.write('{0} null  \n'.format(aln.query_name))


def parse_sam(sam_file):
    """
    Parse the SAM file for alignment lengths, where each
    query can accumulate total length there are multiple
    record against the same reference.

    :param sam_file:
    :return: dictionary of Alignments objects
    """
    align_repo = OrderedDict()
    with open(sam_file, 'r') as h_in:
        for l in h_in:
            # skip header fields
            if l.startswith('@'):
                continue
            fields = l.rsplit()

            aln = Alignment(fields[0], fields[2], count_aligned(fields[5]), seq_length[fields[0]])
            if aln in align_repo:
                align_repo[aln].add_cigar(fields[5])
            else:
                align_repo[aln] = aln
    return align_repo


def parse_last(last_file):
    """
    Parse a LASTAL tabular alignment file for lengths, where
    each alignment of the same query/subject pair accumulate length
    in bases.

    :param last_file: last tabular format alignment file
    :return: dictionary of Alignment objects
    """
    align_repo = OrderedDict()
    with open(last_file, 'r') as h_in:
        for l in h_in:
            # skip header fields
            if l.startswith('#'):
                continue

            fields = l.rsplit()

            qname = fields[6]
            rname = fields[1]
            alen = int(fields[8])
            qlen = int(fields[10])

            # ignore alignment records which fall below mincov with
            # respect to the length of the alignment vs query sequence.
            if float(alen)/float(qlen) < args.mincov:
                continue

            aln = Alignment(qname, rname, alen, qlen)
            if aln in align_repo:
                align_repo[aln].add_bases(alen)
            else:
                align_repo[aln] = aln
    return align_repo


psl_dataline = re.compile(r'^[0-9]+\t')


def parse_psl(psl_file):
    """
    Parse a PSL converted from MAF

    :param psl_file: PSL format alignment file
    :return: dictionary of Alignment objects
    """
    all_hits = 0
    rejected = 0
    align_repo = OrderedDict()
    with open(psl_file, 'r') as h_in:
        for l in h_in:
            all_hits += 1
            # skip header fields
            if not psl_dataline.match(l):
                continue

            all_hits += 1

            fields = l.rsplit()

            qname = fields[9]
            rname = fields[13]
            alen = int(fields[12]) - int(fields[11]) + 1
            qlen = int(fields[10])

            # Taken from BLAT perl script for calculating percentage identity
            matches = int(fields[0])
            mismatches = int(fields[1])
            repmatches = int(fields[2])
            q_num_insert = int(fields[4])
            perid = (1.0 - float(mismatches + q_num_insert) / float(matches + mismatches + repmatches)) * 100.0

            # ignore alignment records which fall below mincov or minid
            # wrt the length of the alignment vs query sequence.
            if float(alen)/float(qlen) < args.mincov or perid < args.minid:
                rejected += 1
                continue

            aln = Alignment(qname, rname, alen, qlen, perid)
            if aln in align_repo:
                align_repo[aln].add_bases(alen)
            else:
                align_repo[aln] = aln

        print 'Rejected {0}/{1} alignments due to constraints on ID {2} and Coverage {3}'.format(
            rejected, all_hits, args.minid, args.mincov)

    return align_repo

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Determine the location of query sequences on the reference')
    parser.add_argument('--simple', default=False, action='store_true',
                        help='Record simple truth table')
    parser.add_argument('--minid', type=float, required=False, default=95.0,
                        help='Minimum percentage identity for alignment')
    parser.add_argument('--minlen', type=int, required=False, default=1000,
                        help='Minimum length in bp')
    parser.add_argument('--mincov', type=float, required=False, default=0.5,
                        help='Minimum coverage of query by alignment')
    parser.add_argument('--afmt', choices=['bwa', 'last', 'psl'], default='last', help='Alignment file format')
    parser.add_argument('--qf', metavar='FASTA', help='Query fasta sequences')
    parser.add_argument('--ofmt', choices=['flat', 'yaml'], default='yaml', help='Output format')
    parser.add_argument('alignment_file', metavar='ALIGNMENT_FILE', nargs=1,
                        help='SAM file or LAST table with query sequences aligned to reference')
    parser.add_argument('output_file', metavar='OUTPUT', nargs=1, help='Output file name')
    args = parser.parse_args()


    align_repo = None
    seq_length = None
    # minimum check that the necessary information has been supplied.
    # if so, fetch lengths of all query sequences
    if args.afmt == 'bwa':
        if args.qf is None:
            raise Exception('BWA alignment format also requires query fasta to be supplied')
        # Calculate the length of query sequences
        seq_length = {rec.id: len(rec) for rec in SeqIO.parse(args.qf, 'fasta')}

        align_repo = parse_sam(args.alignment_file[0])

    elif args.afmt == 'last':
        align_repo = parse_last(args.alignment_file[0])

    elif args.afmt == 'psl':
        align_repo = parse_psl(args.alignment_file[0])

    if args.ofmt == 'flat':
        print 'Soft results always enabled for flat output format'

    print 'Accepted {0} alignments'.format(len(align_repo))

    #
    # We need to decide on assignment.
    #
    # There are query sequences which are assigned to more than one reference.
    # What shall be done with these? Simply taking the largest alignment
    # as the winner might be acceptable for large differences, but I've
    # seen examples where both are significant alignments and picking one
    # is quite arbitrary.
    #

    if not args.simple:
        '''
        Write out all assignments of queries to references.
        '''
        print 'Writing {0} soft (overlapping) assignments of queries to references'.format(len(align_repo.values()))
        write_output(args.simple, args.ofmt, args.minlen, args.mincov,
                     args.minid, align_repo.values(), args.output_file[0])

    else:
        '''
        Pick a winner

        For each scaffold, determine the best alignment by length. The alignment subject
        will then be considered the true source.
        '''
        winners = OrderedDict()
        for aln in align_repo.values():
            if aln.query_name not in winners or winners[aln.query_name].align_length < aln.align_length:
                winners[aln.query_name] = aln

        print 'Reduced to {0} winning hard assignments'.format(len(winners))
        write_output(args.simple, args.ofmt, args.minlen, args.mincov,
                     args.minid, winners.values(), args.output_file[0])
