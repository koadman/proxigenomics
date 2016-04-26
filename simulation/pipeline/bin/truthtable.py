from collections import OrderedDict, Iterable, Counter
import pandas as pd
import copy
import numpy
import yaml

YAML_WIDTH = 1000


# dump OrderedDict like a regular dict. Will not reserialize ordered in this case.
def order_rep(dumper, data):
    return dumper.represent_mapping(u'tag:yaml.org,2002:map', data.items(), flow_style=False)

yaml.add_representer(OrderedDict, order_rep)


class Assignment(yaml.YAMLObject):

    yaml_tag = u'!Assignment'
    yaml_flow_style = True

    def __init__(self, mapping):
        self.mapping = mapping
        self.weight = 1

    def get_classes(self):
        return self.mapping.keys()

    def get_primary_class(self):
        _v = sorted(self.mapping, key=self.mapping.get, reverse=True)
        return _v[0]

    def mean_proportion(self):
        return numpy.mean(self.mapping.values())

    def num_classes(self):
        return len(self.get_classes())

    def __repr__(self):
        return str(self)

    def __str__(self):
        return 'weight={0} mapping={1}'.format(self.weight, str(self.mapping))


class TruthTable(object):
    """
    Class which represents a truth table in the general sense.
    This can be either a prediction or the ground truth.

    The object stores object:class assignments along with a value
    representing the support/strength of that assignment. The intention
    is that this support can be used to convert multi-class assignments
    to a single most significant assignment. It is up to the user
    supply to choose a useful support value.
    """
    def __init__(self):
        self.asgn_dict = {}
        self.label_map = {}
        self.label_count = Counter()

    def __len__(self):
        """
        Length of truth table is equal to the number of assignment classes.
        :return: number of assignment classes
        """
        return len(self.asgn_dict.keys())

    def print_tally(self):
        n_symbol = len(self.label_count)
        n_assignments = sum(self.label_count.values())
        n_objects = len(self.asgn_dict.keys())
        degen_ratio = 100.0 * n_assignments / n_objects - 100.0
        print '{0} symbols in table, {1} assignments of {2} objects ({3:.1f}% degeneracy)'.format(
            n_symbol, n_assignments, n_objects, degen_ratio)
        n_obj = float(sum(self.label_count.values()))
        print 'ext_symb\tint_symb\tcount\tpercentage'
        for ci in sorted(self.label_count, key=self.label_count.get, reverse=True):
            print '{0}\t{1}\t{2}\t{3:5.3f}'.format(ci, self.label_map[ci],
                self.label_count[ci], self.label_count[ci] / n_obj)

    def get_weights(self):
        obj_weights = {}
        for k, v in self.asgn_dict.iteritems():
            obj_weights[k] = v.weight
        return obj_weights

    def object_set(self):
        return set(self.asgn_dict.keys())

    def soft(self, universal=False):
        """
        :return simple dictionary with degenerate classification
        """
        _s = OrderedDict()
        _keys = sorted(self.asgn_dict.keys())
        for k in _keys:
            #_s[k] = sorted(self.asgn_dict[k].get_classes())
            if universal:
                # relabel with universal symbols if requested
                labels = self.asgn_dict[k].get_classes()
                _s[k] = set([self.label_map[l] for l in labels])
            else:
                _s[k] = set(self.asgn_dict[k].get_classes())

        return _s

    def hard(self, universal=False):
        """
        Convert TT to a plain dictionary with the single most significant classification only.
        In the case of a tie, no effort is made to be uniformly random in the case of a tie and
        dependent on the behaviour of sort.
        :return plain dict with only one class->cluster mapping.
        """
        _s = OrderedDict()
        _keys = sorted(self.asgn_dict.keys())
        for _k in _keys:
            #_s[_k] = self.asgn_dict[_k].get_primary_class()
            pc = self.asgn_dict[_k].get_primary_class()
            if universal:
                pc = self.label_map[pc]
            _s[_k] = {pc}
        return _s

    def get(self, key):
        return self.asgn_dict.get(key)

    def put(self, key, value, weight=None):
        self.asgn_dict[key] = Assignment(value)
        if weight:
            self.asgn_dict[key].weight = weight

    def update_from_yaml(self, yd):
        self.asgn_dict = yd
        for asgn_i in yd.values():
            self.label_count.update(asgn_i.mapping.keys())
        labels = sorted(self.label_count.keys())
        self.label_map = dict((l, n) for n, l in enumerate(labels, 1))

    def update(self, dt, weights=None, min_score=0):
        """
        Initialise the assignment dictionary and also generate a mapping of
        class symbol to the positive integers. We can use this as a universal
        symbol basis.

        :param dt: the dictionary to initialise from
        """
        all_asgn = 0
        filt_obj = 0
        filt_asgn = 0
        all_obj = len(dt)

        for k, v in dt.iteritems():
            v_filt = dict((kv, vv) for kv, vv in v.iteritems() if int(vv) >= min_score)
            filt_asgn += len(v) - len(v_filt)
            all_asgn += len(v)
            if len(v_filt) == 0:
                filt_obj += 1
                continue
            self.asgn_dict[str(k)] = Assignment(v_filt)
            if weights:
                self.asgn_dict[str(k)].weight = weights[k]
            self.label_count.update(v_filt.keys())

        if filt_asgn > 0:
            print 'Filtered {0}/{1} assignments and {2}/{3} objects below minimum score {4}'.format(
                filt_asgn, all_asgn, filt_obj, all_obj, min_score)

        labels = sorted(self.label_count.keys())
        self.label_map = dict((l, n) for n, l in enumerate(labels, 1))

    def to_vector(self):
        vec = {}
        keys = sorted(self.asgn_dict.keys())
        hd = self.hard()
        for k in keys:
            vec[k] = hd[k]
        return vec

    @staticmethod
    def _write(dt, pathname, mode='w'):
        """
        Write out bare dictionary.
        :param dt: the dictionary to write
        :param pathname: the path for output file
        """
        with open(pathname, mode) as h_out:
            yaml.dump(dt, h_out, default_flow_style=False, width=YAML_WIDTH)

    def write(self, pathname):
        """
        Write the full table in YAML format"
        :pathname the output path
        """
        TruthTable._write(self.asgn_dict, pathname)

    def write_soft(self, pathname):
        """
        Write a plain dictionary representation in YAML format.
        :pathname the output path
        """
        TruthTable._write(self.soft(), pathname)

    def write_hard(self, pathname):
        """
        Write a plain dictionary representation of only the most significant
        object:class assignments.
        :pathname the output path
        """
        TruthTable._write(self.hard(), pathname)


def read_truth(pathname):
    """
    Read a TruthTable in YAML format
    :param pathname: path to truth table
    :param min_score: ignore objects with scores/weights/lengths below the given threshold.
    :return: truth table
    """
    with open(pathname, 'r') as h_in:
        tt = TruthTable()
        tt.update_from_yaml(yaml.load(h_in))
        return tt


def read_mcl(pathname):
    """
    Read a MCL solution file converting this to a TruthTable
    :param pathname: mcl output file
    :return: truth table
    """

    with open(pathname, 'r') as h_in:
        # read the MCL file, which lists all members of a class on a single line
        # the class ids are implicit, therefore we use line number.
        mcl = {}
        for ci, line in enumerate(h_in, start=1):
            objects = line.rstrip().split()
            for oi in objects:
                if oi not in mcl:
                    mcl[oi] = {}
                mcl[oi][ci] = 1.0  # there are no weights, therefore give them all 1
        # initialise the table
        tt = TruthTable()
        tt.update(mcl)
    return tt


def unique_labels(dt):
    """
    Extract the unique set of class labels used in the truth table
    :param dt: dictionary representation of truth table to analyze
    :return: sorted set of unique labels
    """
    labels = set()
    for v in dt.values():
        if isinstance(v, Iterable):
            labels.update(v)
        else:
            labels.add(v)
    return sorted(labels)


def crosstab(dt1, dt2):
    """
    Cross-tabulate two truth tables on hard clusters.

    :param dt1: first dictionary rep of truth table
    :param dt2: second dictionary rep of truth table
    :return: pandas dataframe
    """
    joined_keys = sorted(set(dt1.keys() + dt2.keys()))

    rows = unique_labels(dt1)
    cols = unique_labels(dt2)
    ctab = pd.DataFrame(0, index=rows, columns=cols)

    for k in joined_keys:
        if k in dt1 and k in dt2:
            i1 = dt1[k]
            i2 = dt2[k]
            ctab.loc[i1, i2] += 1

    return ctab


def simulate_error(tt, p_mut, p_indel, extra_symb=[]):
    """
    Simple method for introducing error in a truth table. This is useful when
    testing clustering metrics (Fm, Vm, Bcubed, etc). By default, the list of possible
    symbols is taken from those already assigned, but a user may provide additional
    symbols. These can provide a useful source of novelty, when for instance
    an object is already assigned to all existing class symbols.

    :param tt: the truth table to add error
    :param p_mut: the probably of a class mutation
    :param p_indel: the probability of deletion or insertion of a class to an object
    :param extra_symb: extra class symbols for inserting
    :return: truth table mutatant
    """
    symbols = list(set(tt.label_count.keys() + extra_symb))
    print symbols
    mut_dict = copy.deepcopy(tt.asgn_dict)

    for o_i in mut_dict.keys():
        others = list(set(symbols) - set(mut_dict[o_i].mapping))
        if numpy.random.uniform() < p_mut:
            if len(others) > 0:
                c_mut = numpy.random.choice(others, 1)[0]
                c_old = numpy.random.choice(mut_dict[o_i].mapping.keys(), 1)[0]

                # retain the weighting from the original to mutated
                # we do this pedantically so its easy to read
                weight = mut_dict[o_i].mapping[c_old]
                mut_dict[o_i].mapping[c_mut] = weight
                del mut_dict[o_i].mapping[c_old]

        if numpy.random.uniform() < p_indel:
            # flip a coin, delete or insert
            if numpy.random.uniform() < 0.5:
                # delete
                c_del = numpy.random.choice(mut_dict[o_i].mapping.keys(), 1)[0]
                del mut_dict[o_i].mapping[c_del]
            elif len(others) > 0:
                # insert from 'others'
                c_add = numpy.random.choice(others, 1)[0]
                num_cl = mut_dict[o_i].num_classes()
                adj_fac = float(num_cl / (num_cl+1.))
                ins_prop = mut_dict[o_i].mean_proportion() * adj_fac
                for k in mut_dict[o_i].mapping:
                    mut_dict[o_i].mapping[k] *= adj_fac
                mut_dict[o_i].mapping[c_add] = ins_prop

    mut_tt = TruthTable()
    mut_tt.update_from_yaml(mut_dict)
    return mut_tt
