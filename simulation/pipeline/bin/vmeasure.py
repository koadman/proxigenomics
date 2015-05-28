#!/usr/bin/env python

from math import log
import truthtable as tt
import sys
import numpy
import pipeline_utils

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
#    print ct_sum
    pi = numpy.sum(ct, axis=1)
    pj = numpy.sum(ct, axis=0)
#    print pi,pj
    outer = numpy.outer(pi, pj)
#    print outer
    nnz = ct != 0.0
#    print nnz
    ct_norm = ct[nnz]
    log_ct = numpy.log(ct[nnz])
    ct_norm /= ct_sum
    #print numpy.log(outer[nnz] / (pi.sum() * pj.sum()))
    log_outer = - numpy.log(outer[nnz]) + log(pi.sum()) + log(pj.sum())
    mi = (ct_norm * (log_ct - log(ct_sum)) + ct_norm * log_outer)
    return mi.sum()


def v_measure(ct):
    H = entropy(ct)
    mutual_info = mutual_information(ct)
    print mutual_info
    homogen = mutual_info / H['C'] if H['C'] else 1.0
    complet = mutual_info / H['K'] if H['K'] else 1.0

    if homogen + complet == 0.0:
        v_measure = 0.0
    else:
        v_measure = (2.0 * homogen * complet) / (homogen + complet)
    return {'homogeneity': float(homogen), 'completeness': float(complet), 'v_measure': float(v_measure)}


'''

The follow functions should replace the existing ones above.

'''

def calculate_entropy(ct):
    N = ct.sum()
    pi = ct.sum(axis=1)
    pj = ct.sum(axis=0)

    s = 0
    for i in xrange(ct.shape[0]):
        if pi[i] > 0.0:
            s += pi[i] / N * log(pi[i] / N)
    HC = - s

    s = 0
    for j in xrange(ct.shape[1]):
        if pj[j] > 0.0:
            s += pj[j] / N * log(pj[j] / N)
    HK = - s

    return HC, HK


def calculate_condentropy(ct):
    N = ct.sum()
    pi = ct.sum(axis=1)
    pj = ct.sum(axis=0)

    s = 0
    for j in xrange(ct.shape[1]):
        if pj[j] > 0.0:
            for i in xrange(ct.shape[0]):
                if ct[i, j] > 0.0:
                    s += ct[i, j] / N * log(ct[i, j] / pj[j])
    HCK = - s

    s = 0
    for i in xrange(ct.shape[0]):
        if pi[i] > 0.0:
            for j in xrange(ct.shape[1]):
                if ct[i, j] > 0.0:
                    s += ct[i, j] / N * log(ct[i, j] / pi[i])
    HKC = - s

    return HCK, HKC

def calculate_score(ct):
    HC, HK = calculate_entropy(ct)
    HCK, HKC = calculate_condentropy(ct)

#    print 'HC={0} HK={1}'.format(HC,HK)
#    print 'HCK={0} HKC={1}'.format(HCK,HKC)

    homog = 1.0 if HCK == 0 else 1.0 - HCK / HC
    compl = 1.0 if HKC == 0 else 1.0 - HKC / HK

    vmeas = 0
    # harmonic mean if denominator is not zero
    if homog + compl != 0:
        vmeas = 2.0 * (homog * compl) / (homog + compl)

    return {'homogeneity': float(homog), 'completeness': float(compl), 'v_measure': float(vmeas)}


if len(sys.argv) != 4:
    print 'Usage [truth] [prediction] [output]'
    sys.exit(1)

truth = tt.read_truth(sys.argv[1])
pred = tt.read_mcl(sys.argv[2])

ct = tt.crosstab(truth.hard(), pred.hard())

# Write the resulting table to stdout
print 'Contigency table [rows=truth, cols=prediction]'
print ct

mat_ct = ct.astype(numpy.float).as_matrix()
# Calculate measures
score = calculate_score(mat_ct)

# Write the scores to the output file
pipeline_utils.write_data(sys.argv[3], score)
