#!/usr/bin/env python

import os
import argparse
import yaml
import pandas
import sys
import traceback

parser = argparse.ArgumentParser(description='Recursively scan path and tabulate scoring data')
parser.add_argument('path', nargs=1, help='Path to scan')
parser.add_argument('csv', nargs=1, help='Output CSV')
args = parser.parse_args()

score_cols = None
lead_cols = None
score_table = []

for path, dirs, files in os.walk(args.path[0]):

    factors = None
    score_vals = {}
    for fn in files:

        if fn.endswith(('.vm', '.f1', '.bc')):

            # variational elements
            factors = path.split('/')

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
                    scores = yaml.load(h_in)
                    if scores is None:
                        sys.stderr.write('{0}/{1} was empty\n'.format(path, fn))
                        continue
                    to_keep = dict(('{0}.{1}'.format(score_method, k), v) for k, v in scores.iteritems())
                    score_vals.update(to_keep)
            except Exception as ex:
                traceback.print_exc(ex)
                raise ex

    if len(score_vals) > 0:
        if score_cols is None:
            score_cols = sorted(score_vals.keys())
        score_table.append(factors + [score_vals.get(sc, 'NA') for sc in score_cols])

df = pandas.DataFrame(score_table, columns=lead_cols + score_cols)
df.to_csv(args.csv[0])
