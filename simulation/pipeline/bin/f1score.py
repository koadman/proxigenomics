#!/usr/bin/env python
from munkres import Munkres, make_cost_matrix
from numpy import diag, average
from numpy.random import randint
from pandas import read_csv, crosstab
from sys import argv, exit, maxsize

def costMatrix(contingencyTable):
    """Hungarian method assumes the goal of minimum cost, while our goal with a
    contigency table is maximum cost. To deal with this, the table is inverted in
    an additive sense.
    """
    return make_cost_matrix(contingencyTable, lambda cost: maxsize - cost)

def matchLabels(contingencyTable):
    """Attempt to match the labels between ground truth and prediction using
    Hungarian algorithm. Extra prediction labels are moved to the right most
    columns.
    """
    m = Munkres()
    ind = m.compute(costMatrix(contingencyTable.as_matrix()))
    return rearrangeColumns(ind, contingencyTable.copy())

def truePositives(contingencyTable):
    """Taken as the diagonal of the matrix"""
    return diag(contingencyTable)

def falseNegatives(contingencyTable):
    """Off diagonal elements in rows"""
    fn = [0] * contingencyTable.shape[0]
    for i in range(contingencyTable.shape[0]):
        for j in range(contingencyTable.shape[1]):
            if i != j:
                fn[i] += contingencyTable.iloc[i,j]
    return fn

def falsePositives(contingencyTable):
    """Off diagonal elements in columns"""
    fp = [0] * contingencyTable.shape[1]
    for j in range(contingencyTable.shape[1]):
        for i in range(contingencyTable.shape[0]):
            if i != j:
                fp[j] += contingencyTable.iloc[i,j]
    return fp

def recall_macro(contingencyTable):
    tp = truePositives(contingencyTable)
    fn = falseNegatives(contingencyTable)
    return average(tp) / (average(tp) + average(fn))

def precision_macro(contingencyTable):
    tp = truePositives(contingencyTable)
    fp = falsePositives(contingencyTable)
    return average(tp) / (average(tp) + average(fp))

def f1_score_macro(contingencyTable):
    """Balanced macro F-measure. This gives equal weight to all classes."""
    tp = truePositives(contingencyTable)
    fn = falseNegatives(contingencyTable)
    fp = falsePositives(contingencyTable)
    return 2.0*average(tp) / (2.0*average(tp) + average(fn) + average(fp))

def rearrangeColumns(indices, contingencyTable):
    """Returns a dataframe where the indices are rearranged to place the maximum value
    on the diagonals. This optimistically matches rows to columns in a classification.
    Additional classes not in the ground truth are left in the right-most columns.
    """
    # sort the table for maximum diagonal elements
    cn_in = contingencyTable.columns.tolist()
    cn_out = [None] * len(cn_in)
    moved = [False] * len(cn_out)
    for i,j in indices:
         cn_out[i] = cn_in[j]
         moved[j] = True
    
    # put losing classes on the end
    last = sum(moved)
    for i,mv in enumerate(moved):
        if not mv:
            cn_out[last] = cn_in[i]
            last += 1
#    print "cout=",cn_out
#    print "cname=",cn_in
#    print "moved=",moved
    return contingencyTable[cn_out]

def overClustered(dataFrame):
    """Check if a matrix has more rows than columns."""
    return dataFrame.shape[0] > dataFrame.shape[1]

def addPaddingColumns(dataFrame):
    """Add dummy columns until matrix is sqaure. We only require this
    for rectangular matrices when there are more rows than columns.
    """
    shp = dataFrame.shape
    nDummy =  shp[0] - shp[1]
    i = 0
    while i < nDummy:
        # pad matrix with zero columns
        dataFrame['dummy' + str(i)] = [0]*shp[0]
        i += 1

if len(argv) != 3:
	print 'Usage: [prediction table] [output]'
	exit(1)

hOut = open(argv[2],'w')

d = read_csv(argv[1],sep=' ')
if len(d) == 0:
	hOut.write('NA NA NA\n')
	exit(0)

ct = crosstab(d['truth'],d['predict'])

if overClustered(ct):
    addPaddingColumns(ct)

mct = matchLabels(ct)

hOut.write("{f1:.4} {recall:.4} {prec:.4}\n".format(
    f1=f1_score_macro(mct),recall=recall_macro(mct),prec=precision_macro(mct)))

hOut.close()
