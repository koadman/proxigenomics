#!/usr/bin/env python
from scipy.stats import hypergeom
import yaml
import numpy
import argparse
import pandas

def intersection(a, i, b, j):
    return numpy.where((a[:, i] * b[:, j] > 0) & (a[:, i] == b[:, j]))[0]


def probability_class_cluster(clazz, klust):
    D = klust.shape[0]
    sum_a = 0
    for i in range(klust.shape[1]):
        for j in range(clazz.shape[1]):
            N = len(intersection(klust, i, clazz, j))
            if N > 0:
                sum_a += hypergeom.pmf(2, D, 2, N)

    sum_b = 0
    for i in range(klust.shape[1]):
        for j in range(clazz.shape[1]):
            for ip in range(i + 1, klust.shape[1]):
                for jp in range(j + 1, clazz.shape[1]):
                    s_ij = intersection(klust, i, clazz, j)
                    s_ipjp = intersection(klust, ip, clazz, jp)
                    N = len(numpy.intersect1d(s_ij, s_ipjp))
                    if N > 0:
                        sum_b += hypergeom.pmf(2, D, 2, N)

    return sum_a - sum_b


def probability_same(assignment):
    D = assignment.shape[0]
    sum_a = 0
    for i in range(assignment.shape[1]):
        N = sum(assignment[:, i])
        if N > 0:
            sum_a += hypergeom.pmf(2, D, 2, N)

    sum_b = 0
    for i in range(assignment.shape[1]):
        for j in range(i + 1, assignment.shape[1]):
            N = len(intersection(assignment, i, assignment, j))
            if N > 0:
                sum_b += hypergeom.pmf(2, D, 2, N)

    return sum_a - sum_b


def gfm_index(clazz, klust):
    """
    Calculate the Generalised Fowlkes-Mallow index for a given clustering of a set classed objects.
    :param clazz: the true class membership (ndarray)
    :param klust: the predicted clustering (ndarray)
    :return: value 0 and 1.
    """


    return probability_class_cluster(clazz, klust) / (probability_same(clazz) * probability_same(klust)) ** 0.5


def assignment_table(soln, labels=None):
    """
    Create a binary assignment table from a given solution (truth or prediction).
    :param soln: solution to a clustering or ground truth
    :return: pandas dataframe assignment table. cows are objects, columns are class/cluster
    """
    if labels is None:
        # first thing, generate the list of labels used in solution
        # if a label was not used, it will not be discovered.
        labels = set()
        for v in soln.values():
            labels |= set(v)
        labels = numpy.array(sorted(labels))

    # tabulate the assignments for each object for each class/cluster
    atab = numpy.zeros((len(soln), len(labels)))
    for i, v in enumerate(soln.values()):
        for vv in v:
            atab[i][numpy.where(labels == vv)[0]] = 1

    # return a pandas table, helpful as it has row and column names.
    return pandas.DataFrame(atab, index=soln.keys(), columns=labels)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Calculate the Generalised Fowlkes-Mallow index")
    parser.add_argument('class_file', metavar='CLASSIFICATION', nargs=1, help='File containing true classification ')
    parser.add_argument('clust_file', metavar='CLUSTERING', nargs=1, help='File containing clustering prediction')
    args = parser.parse_args()

    clazz = None
    with open(args.class_file[0], 'r') as h_in:
        clazz = yaml.load(h_in)

    klust = None
    with open(args.clust_file[0], 'r') as h_in:
        klust = yaml.load(h_in)

    clazz = assignment_table(clazz).as_matrix()
    klust = assignment_table(klust).as_matrix()
    val = gfm_index(clazz, klust)

    print 'gfm {0}'.format(val)
