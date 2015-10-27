#!/usr/bin/env python

import argparse
import sys

parser = argparse.ArgumentParser(description='Parse raw quast information')
parser.add_argument('tsv', metavar='TSV_FILE', help='Input TSV file')

args = parser.parse_args()

with open(args.tsv, 'r') as h_input:

    rows = {}
    cols = set()
    ir = 0

    for n, line in enumerate(h_input):

        if line.startswith('tree'):
            cn = line.strip().split('\t')
            if len(cn) == 0:
                raise IOError('No fields for line {0}'.format(line))

            # header line
            if len(cols) < len(cn):
                for ic in cn:
                    if ic not in cols:
                        cols.add(ic)

            line = h_input.next()
            dn = line.strip().split('\t')
            if len(dn) != len(cn):
                raise IOError('columns do not match for header and data lines. {0}'.format(n))

            # data line
            rows[ir] = dict(zip(cn, dn))
            ir += 1

    sorted_cols = sorted(cols)

    for n, ic in enumerate(sorted_cols):
        sys.stdout.write('{0}'.format(ic))
        if n < len(sorted_cols)-1:
            sys.stdout.write('\t')

    sys.stdout.write('\n')

    for ir in rows.values():
        for n, ic in enumerate(sorted_cols):
            val = ir.get(ic, '')
            sys.stdout.write('{0}'.format(val))
            if n < len(sorted_cols)-1:
                sys.stdout.write('\t')
        sys.stdout.write('\n')
