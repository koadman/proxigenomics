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
    ovl_C = overlap(C[o1], C[o2])
    if ovl_C == 0:
        raise RuntimeError('cardinality of intersection was zero')

    ovl_L = overlap(L[o1], L[o2])

    return float(min(ovl_C, ovl_L)) / float(ovl_C)


def recall(o1, o2, C, L):
    ovl_L = overlap(L[o1], L[o2])
    if ovl_L == 0:
        raise RuntimeError('cardinality of intersection was zero')
    ovl_C = overlap(C[o1], C[o2])
    return float(min(ovl_C, ovl_L)) / float(ovl_L)

def objects_with_overlap(o, C):
    """
    Determine the union of objects which share at least one class (or clustering)
    with an object o. The union also includes o.

    :param o: the subject object
    :param C: the class (or clustering) of objects
    :return: the set of objects overlapping with o
    """
    # class to object map for convenience
    class_to_obj = {}
    for oi, classes in C.iteritems():
        for ci in classes:
            if ci not in class_to_obj:
                class_to_obj[ci] = []
            class_to_obj[ci].append(oi)

    # iterate over o's classes to build the
    # union of all objects sharing a class with o
    D = set([o])  # also include o.
    o_classes = set(C[o])
    for ci in o_classes:
        # objects also belonging to the class
        D |= set(class_to_obj[ci])

    return D

def extended_bcubed_precision(c, g):
    """
    Calculate the overall average extended Bcubed precision for a gold truth and clustering.

    :param c: the gold truth
    :param g: the clustering
    :return: Bcubed_pre
    """
    pre_overall = 0
    shared_obj = intersection(c.keys(), g.keys())  # this might be incorrect.

    if len(shared_obj) <= 0:
        raise RuntimeError('The intersection of truth and prediction was null')

    for obj1 in shared_obj:

        clustering_objects = objects_with_overlap(obj1, g)

        pre = 0
        for obj2 in clustering_objects:
            # only applicable to non-zero intersection of classses
            if obj2 in c and overlap(c[obj1], c[obj2]) > 0:
                pre += multi_precision_or_recall(obj1, obj2, c, g)

        pre_overall += float(pre) / len(clustering_objects)
        #print 'PRE: {:.4f}'.format(float(pre) / len(clustering_objects))

    print '{0} shared objects used in precision'.format(len(shared_obj))
    return pre_overall / len(shared_obj)


def weighted_extended_bcubed_precision(w, c, g):
    """
    Calculate the overall average extended Bcubed precision for a gold truth and clustering.

    :param c: the gold truth
    :param g: the clustering
    :return: Bcubed_pre
    """
    pre_overall = 0
    shared_obj = intersection(c.keys(), g.keys())  # this might be incorrect.

    if len(shared_obj) <= 0:
        raise RuntimeError('The intersection of truth and prediction was null')

    for obj1 in shared_obj:

        clustering_objects = objects_with_overlap(obj1, g)

        pre = 0
        for obj2 in clustering_objects:
            # only applicable to non-zero intersection of classses

            if obj2 in c:
                ovl_c = overlap(c[obj1], c[obj2])
                if ovl_c > 0:
                    ovl_g = overlap(g[obj1], g[obj2])
                    pre += w[obj2] * float(min(ovl_c, ovl_g)) / float(ovl_c)

        weight_clust = float(sum([w[k] for k in clustering_objects]))
        pre_overall += float(pre) / weight_clust
        #print 'PRE: {:.4f}'.format(float(pre) / weight_clust)

    print '{0} shared objects used in precision'.format(len(shared_obj))
    return pre_overall / len(shared_obj)


def extended_bcubed_recall(c, g):
    """
    Calculate the overall average extended Bcubed recall for a gold truth and clustering.

    :param c: the gold truth
    :param g: the clustering
    :return: Bcubed_rec
    """
    rec_overall = 0
    shared_obj = intersection(c.keys(), g.keys())

    if len(shared_obj) <= 0:
        raise RuntimeError('The intersection of truth and prediction was null')

    for obj1 in shared_obj:

        class_objects = objects_with_overlap(obj1, c)

        rec = 0
        for obj2 in class_objects:
            # only applicable to non-zero intersection of classses
            if obj2 in g and overlap(g[obj1], g[obj2]) > 0:
                rec += multi_precision_or_recall(obj1, obj2, g, c)

        rec_overall += float(rec) / len(class_objects)
        #print 'REC: {:.4f}'.format(float(rec) / len(class_objects))

    print '{0} shared objects used in recall'.format(len(shared_obj))
    return rec_overall / len(shared_obj)


def weighted_extended_bcubed_recall(w, c, g):
    """
    Calculate the overall average extended Bcubed recall for a gold truth and clustering.

    :param c: the gold truth
    :param g: the clustering
    :return: Bcubed_rec
    """
    rec_overall = 0
    shared_obj = intersection(c.keys(), g.keys())

    if len(shared_obj) <= 0:
        raise RuntimeError('The intersection of truth and prediction was null')

    for obj1 in shared_obj:

        class_objects = objects_with_overlap(obj1, c)

        rec = 0
        for obj2 in class_objects:
            # only applicable to non-zero intersection of classses
            if obj2 in g:
                ovl_g = overlap(g[obj1], g[obj2])
                if ovl_g > 0:
                    ovl_c = overlap(c[obj1], c[obj2])
                    rec += w[obj2] * float(min(ovl_c, ovl_g)) / float(ovl_g)

        weight_class = float(sum([w[k] for k in class_objects]))
        rec_overall += float(rec) / weight_class
        #print 'REC: {:.4f}'.format(float(rec) / weight_class)

    print '{0} shared objects used in recall'.format(len(shared_obj))
    return rec_overall / len(shared_obj)


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


def weighted_extended_bcubed(w, c, g):
    """
    Calculate the overall average extended Bcubed F measure for a gold truth and clustering.
    The F measure is the harmonic mean of the extended Bcubed precision and recall.

    :param c: the gold truth
    :param g: the clustering
    :return: precision, recall and F.
    """
    pre = weighted_extended_bcubed_precision(w, c, g)
    rec = weighted_extended_bcubed_recall(w, c, g)
    fbcubed = 2.0 * pre * rec / (pre + rec)
    return {'pre': pre, 'rec': rec, 'f': fbcubed}


def write_msg(filename, msg):
    if filename is None:
        print msg
    else:
        with open(filename, 'w') as hout:
            hout.write(msg + '\n')


if __name__ == '__main__':

    def count_labels(cl):
        cl_names = Counter()
        for v in cl.values():
            cl_names.update(v)
        return cl_names

    from collections import Counter
    import truthtable as tt
    import pipeline_utils
    import argparse
    import sys

    parser = argparse.ArgumentParser(description='Calculate extended bcubed metric')
    parser.add_argument('--drop-small', dest='small_thres', metavar='THRES_PROPORTION', type=float, help='Filter clusters which fall belong minimum proportion ')
    #parser.add_argument('--weighted', dest='weight_csv', metavar='WEIGHT_CSV', help='Calculate a weighted bcubed')
    parser.add_argument('--weighted', default=False, action='store_true', help='Use weighted terms')
    parser.add_argument('-o', '--output', help='Output file')
    parser.add_argument('truth', metavar='TRUTH', nargs=1, help='Truth table (yaml format)')
    parser.add_argument('pred', metavar='PREDICTION', nargs=1, help='Prediction table (mcl format)')
    args = parser.parse_args()

    try:
        # read truth and convert to basic soft table
        truth = tt.read_truth(args.truth[0])
        if len(truth) == 0:
            raise RuntimeWarning('Truth table contains no assignments: {0}'.format(args.truth[0]))

        print 'Truth Statistics'
        truth.print_tally()
        if args.weighted:
            weights = truth.weights()

        truth = truth.soft(True)

        # read clustering and convert to basic soft table
        clustering = tt.read_mcl(args.pred[0])
        if len(clustering) == 0:
            raise RuntimeWarning('Clustering contains no assignments: {0}'.format(args.pred[0]))
        if args.small_thres:
            clustering.filter_class(args.small_thres)

        print 'Clustering Statistics'
        clustering.print_tally()
        clustering = clustering.soft(True)

    except RuntimeWarning as wn:
        write_msg(args.output, wn.message)
        sys.exit(0)

    # initialise the output stream
    if args.output is None:
        args.output = sys.stdout
    else:
        args.output = open(args.output, 'w')

    # Read weights from CSV.
    #weights = None
    #if args.weight_csv is not None:
    #    weights = {}
    #    with open(args.weight_csv, 'r') as h_in:
    #        for line in h_in:
    #            fields = line.strip().split()
    #            if len(fields) != 2:
    #                raise IOError('weight csv did not contain 2 columns')
    #            if fields[0] in weights:
    #                raise IOError('weight csv contains duplicate keys')
    #            try:
    #                weights[fields[0]] = float(fields[1])
    #            except ValueError as er:
    #                sys.stderr.write('Warning: caught conversion error for node weights. [{0}]\n'.format(er.message))
    if args.weighted:
        pipeline_utils.write_to_stream(args.output, weighted_extended_bcubed(weights, truth, clustering))

    else:
    #unweighted bcubed
        pipeline_utils.write_to_stream(args.output, extended_bcubed(truth, clustering))

