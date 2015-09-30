#!/usr/bin/env python

import os
import argparse
import yaml
import pandas
import sys
import traceback


def read_yaml(h_in):
    return yaml.load(h_in)


def read_adhoc(h_in):
    scores = {}
    for line in h_in:
        if line.startswith('Eigen value'):
            scores['val'] = line.split()[3]
    return scores


parser = argparse.ArgumentParser(description='Recursively scan path and tabulate scoring data')
parser.add_argument('--suffix', nargs='*', default=['.vm', '.f1', '.bc'],
                    help='One or more filename suffixes to target in recursive search [(.vm, .f1, .bc)]')
parser.add_argument('-y', '--yaml-vals', action="store_true", default=False,
                    help='Scores are stored in YAML format files')
parser.add_argument('-p', '--path', metavar='DIR', required=True, help='Input path to scan')
parser.add_argument('-o', '--output', metavar='FILE', required=True, help='Output table')
args = parser.parse_args()

print args
score_cols = None
lead_cols = None
score_table = []

suffixes = tuple(args.suffix)

if args.yaml_vals:
    read_vals = read_yaml
else:
    read_vals = read_adhoc

for path, dirs, files in os.walk(args.path):

    factors = None
    score_vals = {}
    for fn in files:

        if fn.endswith(suffixes):

            # variational elements
            factors = path.split(os.path.sep)

            # leading columns are just indexed
            if lead_cols is None:
                lead_cols = ['c{}'.format(n) for n, v in enumerate(factors, start=1)]

            elif len(lead_cols) != len(factors):
                raise RuntimeError('There are different variational levels within this path hierarchy\n' +
                '{0}\n{1}\n{2}'.format(path, lead_cols, factors))

            # scores
            try:
                with open(os.path.join(path, fn), 'r') as h_in:
                    score_method = os.path.splitext(fn)[1][1:]
                    scores = read_vals(h_in)
                    if scores is None:
                        sys.stderr.write('{0}/{1} was empty\n'.format(path, fn))
                        continue
                    to_keep = dict(('{0}.{1}'.format(score_method, k), v) for k, v in scores.iteritems())
                    score_vals.update(to_keep)
            except Exception as ex:
                print '{0}/{1}'.format(path, fn)
                traceback.print_exc(ex)
                raise ex

    if len(score_vals) > 0:
        if score_cols is None:
            score_cols = sorted(score_vals.keys())
        score_table.append(factors + [score_vals.get(sc, 'NA') for sc in score_cols])

df = pandas.DataFrame(score_table, columns=lead_cols + score_cols)
df.to_csv(args.output)
