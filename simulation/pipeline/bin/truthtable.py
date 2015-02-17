import pandas as pd
import yaml


class TruthTable(dict):
    """
    Class which represents a truth table in the general sense.
    This can be either a prediction or the ground truth.

    The object stores object:class assignments along with a value
    representing the support/strength of that assignment. The intention
    is that this support can be used to convert multi-class assignments
    to a single most significant assignment. It is up to the user
    supply to choose a useful support value.
    """
    def __init__(self, *args, **kwargs):
        super(TruthTable, self).__init__(*args, **kwargs)

    def soft(self):
        """
        :return plain dictionary with degenerate classification
        """
        _s = {}
        _keys = sorted(self.keys())
        for k in _keys:
            _s[k] = sorted(self[k].keys())
        return _s

    def hard(self):
        """
        :return plain dictionary with the single most significant classification only.
        """
        _s = {}
        _keys = sorted(self.keys())
        for _k in _keys:
            _v = self[_k]
            _v = sorted(_v, key=_v.get, reverse=True)
            _s[_k] = _v[0]
        return _s

    @staticmethod
    def _write(table, pathname):
        with open(pathname, 'w') as h_out:
            yaml.dump(table, h_out)

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
        return yaml.load(h_in)


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
                mcl[o] = {n: None}
    return mcl


def crosstab(tt1, tt2):
    """
    Cross-tabulate two truth tables on hard clusters. Call tt.hard() first.

    :param tt1: first truth table
    :param tt2: second truth table
    :return: pandas dataframe
    """
    joined_keys = sorted(set(tt1.keys() + tt2.keys()))
    rows = sorted(set(tt1.values()))
    cols = sorted(set(tt2.values()))
    ctab = pd.DataFrame(0, index=rows, columns=cols)

    for k in joined_keys:
        if k in tt1 and k in tt2:
            i1 = tt1[k]
            i2 = tt2[k]
            ctab.loc[i1, i2] += 1

    return ctab


def ordereddict_rep(dumper, data):
    """
    Helper method to handling OrderedDict types in YAML dump.

    :param dumper:
    :param data:
    :return:
    """
    return dumper.represent_mapping(u'tag:yaml.org,2002:map', data.items(), flow_style=False)

