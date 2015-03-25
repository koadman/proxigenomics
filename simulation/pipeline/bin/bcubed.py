#!/usr/bin/env python
# coding=utf-8


"""
Calculate the intersection of two sets, first placing the iterable into a set.
"""
intersection = lambda set_a, set_b: set(set_a) & set(set_b)

"""
Calculate the overlap between two sets. That is, the cardinality of the intersection of set A and set B.
"""
overlap = lambda set_a, set_b: len(intersection(set_a, set_b))

def multi_precision_or_recall(o1, o2, C, L):
    """
    Calculate multiplicity B-cubed precision or recall. These two equations differ only
    in their normalisation, where precision depends on C, recall depends on L. Therefore
    the method emits either by switching C and L at call time.

    Supplying gold truth to C results in precision, while supplying clustering to C results in recall.

    Using notation from:
    Perez-Suarez, A.; Martinez-Trinidad, J. F.; Carrasco-Ochoa, J. A.; Medina-Pagola, J. E. OClustR:
    A new graph-based algorithm for overlapping clustering. Neurocomputing 2013, 121, 234–247.


    :param o1: an object to compare overlap with another
    :param o2: another object
    :param C: the gold truth (precision), the clustering (recall)
    :param L: the clustering (precision), the gold truth (recall)
    :return:
    """
    ovl_C = overlap(C[o1], C[o2]) #len(set(C[o1]) & set(C[o2]))
    if ovl_C == 0:
        raise RuntimeError('cardinality of intersection was zero')

    ovl_L = overlap(L[o1], L[o2]) #len(set(L[o1]) & set(L[o2]))

    return float(min(ovl_C, ovl_L)) / float(ovl_C)


def extended_bcubed_precision(c, g):
    """
    Calculate the overall average extended Bcubed precision for a gold truth and clustering.

    :param c: the gold truth
    :param g: the clustering
    :return: Bcubed_pre
    """
    pre_overall = 0
    n_overall = 0
    shared_obj = intersection(c.keys(), g.keys())  #set(c.keys()) & set(g.keys())
    for obj1 in shared_obj:
        n = 0
        pre = 0
        for obj2 in shared_obj:
            # count the number of objects which share at least one cluster with obj1
            n += overlap(g[obj1], g[obj2])
            # only applicable to non-zero intersection of classses
            if obj2 in g and overlap(c[obj1], c[obj2]) > 0:
                pre += multi_precision_or_recall(obj1, obj2, c, g)
        pre_overall += float(pre)/n
        n_overall += 1
    print 'n_overall={0}'.format(n_overall)
    return pre_overall/n_overall


def extended_bcubed_recall(c, g):
    """
    Calculate the overall average extended Bcubed recall for a gold truth and clustering.

    :param c: the gold truth
    :param g: the clustering
    :return: Bcubed_rec
    """
    rec_overall = 0
    n_overall = 0
    shared_obj = intersection(c.keys(), g.keys())  #set(c.keys()) & set(g.keys())
    for obj1 in shared_obj:
        n = 0
        rec = 0
        for obj2 in shared_obj:
            # count the number of objects which share at least one class with obj1
            n += overlap(c[obj1], c[obj2])
            # only applicable to non-zero intersection of classses
            if overlap(g[obj1], g[obj2]) > 0:
                rec += multi_precision_or_recall(obj1, obj2, g, c)
                #n += 1
        rec_overall += float(rec)/n
        n_overall += 1
    print 'n_overall={0}'.format(n_overall)
    return rec_overall/n_overall


def extended_bcubed(c, g):
    """
    Calculate the overall average extended Bcubed F measure for a gold truth and clustering.
    The F measure is the harmonic mean of the extended Bcubed precision and recall.

    :param c: the gold truth
    :param g: the clustering
    :return: precision, recall and F.
    """
    pre = extended_bcubed_precision(c, g)
    rec = extended_bcubed_recall(c, g)
    fbcubed = 2.0 * pre * rec / (pre + rec)
    return {'pre': pre, 'rec': rec, 'f': fbcubed}

if __name__ == '__main__':

    def count_labels(cl):
        cl_names = Counter()
        for v in cl.values():
            cl_names.update(v)
        return cl_names

    from collections import Counter
    import truthtable as tt
    import pipeline_utils
    import sys

    if len(sys.argv) != 4:
        print 'Usage [truth] [prediction] [output]'
        sys.exit(1)

    # read truth and convert to basic soft table
    truth = tt.read_truth(sys.argv[1])
    print 'Truth Statistics'
    truth.print_tally()
    truth = truth.soft(True)

    # read clustering and convert to basic soft table
    clustering = tt.read_mcl(sys.argv[2])
    print 'Clustering Statistics'
    clustering.print_tally()
    clustering = clustering.soft(True)

    pipeline_utils.write_data(sys.argv[3], extended_bcubed(truth, clustering))
