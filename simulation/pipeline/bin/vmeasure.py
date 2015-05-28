#!/usr/bin/env python

from math import log
import truthtable as tt
import sys
import numpy
import pipeline_utils
import argparse

def entropy(ct):
    """Calculate the maximal entropy of classes C and clusterings K from
    a contingency table.
    Returns dict {C, K}
    """
    ki = ct.sum(axis=0)
    ki_sum = numpy.sum(ki)
    HK = -numpy.sum((ki / ki_sum) * (numpy.log(ki) - log(ki_sum)))
    ci = ct.sum(axis=1)
    ci_sum = numpy.sum(ci)
    HC = -numpy.sum((ci / ci_sum) * (numpy.log(ci) - log(ci_sum)))
    return {'C': HC, 'K': HK}


def mutual_information(ct):
    """Calculate the mutual information score of a contingency table of truth vs prediction
    """
    ct = ct.astype(numpy.float).as_matrix()
    ct_sum = numpy.sum(ct)
    pi = numpy.sum(ct, axis=1)
    pj = numpy.sum(ct, axis=0)
    outer = numpy.outer(pi, pj)
    nnz = ct != 0.0
    ct_norm = ct[nnz]
    log_ct = numpy.log(ct[nnz])
    ct_norm /= ct_sum
    log_outer = - numpy.log(outer[nnz]) + log(pi.sum()) + + log(pj.sum())
    mi = (ct_norm * (log_ct - log(ct_sum)) + ct_norm * log_outer)
    return mi.sum()


def v_measure(ct):
    H = entropy(ct)
    mutual_info = mutual_information(ct)

    homogen = mutual_info / H['C'] if H['C'] else 1.0
    complet = mutual_info / H['K'] if H['K'] else 1.0

    if homogen + complet == 0.0:
        v_measure = 0.0
    else:
        v_measure = (2.0 * homogen * complet) / (homogen + complet)
    return {'homogeneity': float(homogen), 'completeness': float(complet), 'v_measure': float(v_measure)}


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Calculate V_measure')
    parser.add_argument('truth', metavar='TRUTH', nargs=1, help='Truth table (yaml format)')
    parser.add_argument('pred', metavar='PREDICTION', nargs=1, help='Prediction table (mcl format)')
    parser.add_argument('output', metavar='OUTPUT', nargs='?',
                        type=argparse.FileType('w'), default=sys.stdout, help='Output file')
    args = parser.parse_args()

    truth = tt.read_truth(args.truth[0])
    pred = tt.read_mcl(args.pred[0])

    ct = tt.crosstab(truth.hard(), pred.hard())

    # Write the resulting table to stdout
    print 'Contigency table [rows=truth, cols=prediction]'
    print ct

    # Calculate measures
    score = v_measure(ct)

    # Write the scores to the output file
    pipeline_utils.write_to_stream(args.output, score)
