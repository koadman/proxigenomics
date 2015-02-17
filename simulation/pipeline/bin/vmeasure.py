#!/usr/bin/env python

from math import log
import truthtable as tt
import sys
import numpy


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
    return {'homogeneity': homogen, 'completeness': complet, 'v_measure': v_measure}


if len(sys.argv) != 4:
    print 'Usage [truth] [prediction] [output]'
    sys.exit(1)

truth = tt.read_truth(sys.argv[1])
pred = tt.read_mcl(sys.argv[2])

ct = tt.crosstab(truth.hard(), pred.hard())

print 'Contigency table [rows=truth, cols=prediction]'
print ct

with open(sys.argv[3], 'w') as h_out:
    h_out.write('{0[homogeneity]:.4} {0[completeness]:.4} {0[v_measure]:.4}\n'.format(v_measure(ct)))


