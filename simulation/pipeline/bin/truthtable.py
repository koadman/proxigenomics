import pandas as pd
from collections import OrderedDict, Iterable
import yaml


class TruthTable:
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

    def soft(self):
        """
        :return plain dictionary with degenerate classification
        """
        _s = OrderedDict()
        _keys = sorted(self.asgn_dict.keys())
        for k in _keys:
            _s[k] = sorted(self.asgn_dict[k].keys())
        return _s

    def hard(self):
        """
        :return plain dictionary with the single most significant classification only.
        """
        _s = OrderedDict()
        _keys = sorted(self.asgn_dict.keys())
        for _k in _keys:
            _v = self.asgn_dict[_k]
            _v = sorted(_v, key=_v.get, reverse=True)
            _s[_k] = _v[0]
        return _s

    def get(self, key):
        return self.asgn_dict.get(key)

    def put(self, key, value):
        self.asgn_dict[key] = value

    def update(self, dt):
        for k, v in dt.iteritems():
            self.asgn_dict[k] = v

    def add(self, key, value):
        self.asgn_dict[key] = value

    def to_vector(self):
        vec = {}
        keys = sorted(self.asgn_dict.keys())
        hd = self.hard()
        for k in keys:
            vec[k] = hd[k]
        return vec


    @staticmethod
    def _write(table, pathname):
        with open(pathname, 'w') as h_out:
            yaml.dump(table.asgn_dict, h_out, default_flow_style=False)

    def write(self, pathname):
        """
        Write the full table in YAML format"
        :pathname the output path
        """
        TruthTable._write(self, pathname)

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
    :return: truth table
    """
    with open(pathname, 'r') as h_in:
        tt = TruthTable()
        tt.update(yaml.load(h_in))
        return tt


def read_mcl(pathname):
    """
    Read a MCL solution file converting this to a TruthTable
    :param pathname: mcl output file
    :return: truth table
    """
    with open(pathname, 'r') as h_in:
        mcl = TruthTable()
        for n, line in enumerate(h_in, start=1):
            objs = line.rstrip().split()
            for o in objs:
                mcl.add(o, {n: None})
    return mcl


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
