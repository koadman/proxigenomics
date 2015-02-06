#!/usr/bin/env python

from munkres import Munkres, make_cost_matrix
from numpy import diag, average
from pandas import read_csv, crosstab
import sys


def cost_matrix(contingency_table):
    """Hungarian method assumes the goal of minimum cost, while our goal with a
    contigency table is maximum cost. To deal with this, the table is inverted in
    an additive sense.
    """
    return make_cost_matrix(contingency_table, lambda cost: sys.maxsize - cost)


def match_labels(contingency_table):
    """Attempt to match the labels between ground truth and prediction using
    Hungarian algorithm. Extra prediction labels are moved to the right most
    columns.
    """
    m = Munkres()
    ind = m.compute(cost_matrix(contingency_table.as_matrix()))
    return rearrange_columns(ind, contingency_table.copy())


def true_positives(contingency_table):
    """Taken as the diagonal of the matrix"""
    return diag(contingency_table)


def false_negatives(contingency_table):
    """Off diagonal elements in rows"""
    fn = [0] * contingency_table.shape[0]
    for i in range(contingency_table.shape[0]):
        for j in range(contingency_table.shape[1]):
            if i != j:
                fn[i] += contingency_table.iloc[i, j]
    return fn


def false_positives(contingency_table):
    """Off diagonal elements in columns"""
    fp = [0] * contingency_table.shape[1]
    for j in range(contingency_table.shape[1]):
        for i in range(contingency_table.shape[0]):
            if i != j:
                fp[j] += contingency_table.iloc[i, j]
    return fp


def recall_macro(contingency_table):
    tp = true_positives(contingency_table)
    fn = false_negatives(contingency_table)
    return average(tp) / (average(tp) + average(fn))


def precision_macro(contingency_table):
    tp = true_positives(contingency_table)
    fp = false_positives(contingency_table)
    return average(tp) / (average(tp) + average(fp))


def f1_score_macro(contingency_table):
    """Balanced macro F-measure. This gives equal weight to all classes."""
    tp = true_positives(contingency_table)
    fn = false_negatives(contingency_table)
    fp = false_positives(contingency_table)
    return 2.0 * average(tp) / (2.0 * average(tp) + average(fn) + average(fp))


def rearrange_columns(indices, contingency_table):
    """Returns a dataframe where the indices are rearranged to place the maximum value
    on the diagonals. This optimistically matches rows to columns in a classification.
    Additional classes not in the ground truth are left in the right-most columns.
    """
    # sort the table for maximum diagonal elements
    cn_in = contingency_table.columns.tolist()
    cn_out = [None] * len(cn_in)
    moved = [False] * len(cn_out)
    for i, j in indices:
        cn_out[i] = cn_in[j]
        moved[j] = True

    # put losing classes on the end
    last = sum(moved)
    for i, mv in enumerate(moved):
        if not mv:
            cn_out[last] = cn_in[i]
            last += 1
            # print "cout=",cn_out
        #    print "cname=",cn_in
        #    print "moved=",moved
    return contingency_table[cn_out]


def over_clustered(dataFrame):
    """Check if a matrix has more rows than columns."""
    return dataFrame.shape[0] > dataFrame.shape[1]


def add_padding_columns(dataFrame):
    """Add dummy columns until matrix is sqaure. We only require this
    for rectangular matrices when there are more rows than columns.
    """
    shp = dataFrame.shape
    nDummy = shp[0] - shp[1]
    i = 0
    while i < nDummy:
        # pad matrix with zero columns
        dataFrame['dummy' + str(i)] = [0] * shp[0]
        i += 1


if len(sys.argv) != 3:
    print 'Usage: [prediction table] [output]'
    exit(1)

with open(sys.argv[2], 'w') as h_out:

    d = read_csv(sys.argv[1], sep=' ')
    if len(d) == 0:
        h_out.write('NA NA NA\n')
        sys.exit(0)

    ct = crosstab(d['truth'], d['predict'])

    if over_clustered(ct):
        add_padding_columns(ct)

    mct = match_labels(ct)

    h_out.write("{f1:.4} {recall:.4} {prec:.4}\n".format(
        f1=f1_score_macro(mct), recall=recall_macro(mct), prec=precision_macro(mct)))
