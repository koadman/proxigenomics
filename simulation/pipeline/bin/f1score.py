#!/usr/bin/env python

from munkres import Munkres, make_cost_matrix
import pandas as pd
import numpy as np
import truthtable as tt
import pipeline_utils
import argparse
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
    return np.diag(contingency_table)


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
    return float(np.average(tp) / (np.average(tp) + np.average(fn)))


def precision_macro(contingency_table):
    tp = true_positives(contingency_table)
    fp = false_positives(contingency_table)
    return float(np.average(tp) / (np.average(tp) + np.average(fp)))


def f1_score_macro(contingency_table):
    """Balanced macro F-measure. This gives equal weight to all classes."""
    tp = true_positives(contingency_table)
    fn = false_negatives(contingency_table)
    fp = false_positives(contingency_table)
    return float(2.0 * np.average(tp) / (2.0 * np.average(tp) + np.average(fn) + np.average(fp)))


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


def print_table(df):
    """
    Create a temporary version of the datatable which includes marginal sums of
    columns and rows. Uses the original column and row names.
    :param df: pandas dataframe to print
    """
    a = np.empty((df.shape[0], df.shape[1]+1), dtype=int)
    a[:, 1:] = df
    a[:, 0] = np.sum(df, axis=1)
    b = np.empty((a.shape[0]+1, a.shape[1]), dtype=int)
    b[1:, :] = a
    b[0, :] = np.sum(a, axis=0)
    print pd.DataFrame(b,
                       columns=['Sum'] + df.columns.values.tolist(),
                       index=['Sum'] + df.index.values.tolist())


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Calculate F1 metric')
    parser.add_argument('-s', '--min-score', type=int, help='minimum truth object score', default=0)
    parser.add_argument('truth', nargs=1, help='truth table in yaml format')
    parser.add_argument('pred', nargs=1, help='prediction in MCL format')
    parser.add_argument('output', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='Output file')
    args = parser.parse_args()

    print 'Reading truth table...'
    truth = tt.read_truth(args.truth[0], args.min_score)
    print 'Reading prediction...'
    pred = tt.read_mcl(args.pred[0])

    print 'Creating contingency table...'
    ct = tt.crosstab(truth.hard(), pred.hard())

    print
    print 'Contigency table [rows=truth, cols=prediction] contains {0} elements'.format(ct.shape[0] * ct.shape[1])
    print_table(ct)
    print

    if over_clustered(ct):
        add_padding_columns(ct)
        print 'Squaring table with dummy classes'
        print_table(ct)
        print

    # Write the table to stdout
    print 'Matching labels using Munkres, algorithm suffers from high polynomial order...'
    mct = match_labels(ct)
    print 'Aligned contigency table'
    print_table(mct)

    # Calculate measure
    score = {'f1': f1_score_macro(mct),
             'recall': recall_macro(mct),
             'prec': precision_macro(mct)}

    # Write measure to file
    pipeline_utils.write_to_stream(args.output, score)
