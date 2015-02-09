#!/usr/bin/env python
from scipy.stats import hypergeom
import numpy


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
                    Sij = intersection(klust, i, clazz, j)
                    Sipjp = intersection(klust, ip, clazz, jp)
                    N = len(numpy.intersect1d(Sij, Sipjp))
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
    :param clazz: the true class membership
    :param klust: the predicted clustering
    :return: value 0 and 1.
    """
    return probability_class_cluster(clazz, klust) / (probability_same(clazz) * probability_same(klust)) ** 0.5


c = numpy.array(
    [[1, 0, 1], [1, 0, 0], [1, 0, 0], [0, 1, 1], [0, 1, 0],
     [1, 0, 0], [0, 1, 0], [0, 1, 1], [0, 1, 1], [0, 1, 1]])

k = numpy.array(
    [[1, 0], [1, 0], [1, 0], [1, 0], [1, 0],
     [0, 1], [0, 1], [0, 1], [0, 1], [0, 1]])


print gfm_index(c, k)
