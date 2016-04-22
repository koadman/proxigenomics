#!/usr/bin/env python
import sys

if len(sys.argv) != 2:
    print 'usage: [stat file]'
    sys.exit(1)

md = []
seq_sum = 0
len_sum = 0
n_line = 0

with open(sys.argv[1], 'r') as input_h:

    for l in input_h:
        l = l.strip()
        if not l:
            continue

        tok = l.split()

        if l.startswith('##'):
            md.append(tok[1])
        elif l.startswith('N'):
            continue
        else:
            if len(tok) != 6:
                sys.stderr.write('{0} assuming end of tablular data after {1} rows\n'.format(sys.argv[1], n_line))
                break
            n_line += 1
            seq_sum += int(tok[2])
            len_sum += int(tok[3])
                

print sys.argv[1],' '.join(md),seq_sum,len_sum    
