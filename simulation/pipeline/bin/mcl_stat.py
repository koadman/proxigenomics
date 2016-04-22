#!/usr/bin/env python
import truthtable as tt
import networkx as nx
import argparse
import sys

def tally_truth(t, lengths):
    tally = {ti: {'seq_num':0, 'length':0, 'max_len':0} for ti in t.label_count}
    t = t.soft()    
    for ni, asigns in t.iteritems():
        for ai in asigns:
            try:
                tally[ai]['length'] += lengths[ni]
                tally[ai]['seq_num'] += 1
                if lengths[ni] > tally[ai]['max_len']:
                    tally[ai]['max_len'] = lengths[ni]
            except BaseException:
                import traceback
                exc_type, exc_value, exc_traceback = sys.exc_info()                print "*** print_exception:"                traceback.print_exception(exc_type, exc_value, exc_traceback, file=sys.stderr)
                sys.stderr.write('args: {0}\n'.format(args))
                sys.exit(1)

    return tally

parser = argparse.ArgumentParser(description='Calculate statistics on clustering result')
parser.add_argument('--min-prop', default='0', type=float, help='Exclude groupings which contain less than this proportion of objects [0]')
parser.add_argument('--cfmt', default='mcl', choices=['mcl','yaml'], help='Input clustering file format [mcl]')
parser.add_argument('--lfmt', default='graphml', choices=['graphml','idxstats'], help='Input lengths in graph or bam idxstats format')
parser.add_argument('graph', metavar='LENGTH_DAT', help='Input file in either MCL or YAML format')
parser.add_argument('clustering', metavar='CLUSTERING', help='Input file in either MCL or YAML format')
args = parser.parse_args()

lengths = {}
if args.lfmt == 'graphml':
    g = nx.read_graphml(args.graph)
    for ni in g.nodes_iter():
        lengths[ni] = g.node[ni]['length']
    del g
elif args.lfmt == 'idxstats':
    with open(args.graph, 'r') as in_h:
        n_line = 0
        while True:
            line = in_h.readline()
            n_line += 1
            if line.startswith('*'):
                break
            if not line:
                continue
            tok = line.split()
            if len(tok) != 4:
                raise IOError('bad idxstats file; wrong number of columns at line: {0}'.format(n_line))
            try:
                lengths[tok[0]] = int(tok[1])
            except BaseException as e:
                import traceback
                exc_type, exc_value, exc_traceback = sys.exc_info()                print "*** print_exception:"                traceback.print_exception(exc_type, exc_value, exc_traceback, limit=2, file=sys.stderr)
                sys.stderr.write('args: {0}\n'.format(args))
                sys.exit(1)


if args.cfmt == 'mcl':
    t = tt.read_mcl(args.clustering)
elif args.cfmt == 'yaml':
    t = tt.read_truth(args.clustering)

if t.num_objects() == 0:
    print 'No objects found in {0}'.format(args.clustering)
    sys.exit(1)

try:
    t.filter_extent(args.min_prop, lengths)
except BaseException as e:
    import traceback
    exc_type, exc_value, exc_traceback = sys.exc_info()    print "*** print_exception:"    traceback.print_exception(exc_type, exc_value, exc_traceback, limit=2, file=sys.stderr)
    sys.stderr.write('args: {0}\n'.format(args))

#t.filter_class(args.min_prop)
print '##minimum_proportion {0}'.format(args.min_prop)
print '##global_degeneracy {0}'.format(t.degeneracy(lengths))
print '##mean_overlap {0}'.format(t.mean_overlap(lengths))

tally = tally_truth(t, lengths)

idx = sorted(tally, key=lambda x: (tally[x]['length'],-tally[x]['seq_num']), reverse=True)

print 'N\tCLID\tSEQ_NUM\tLENGTH\tMAX_LEN\tMEAN_LEN'
for n, clid in enumerate(idx,start=1):
    print '{0}\t{1}\t{2}\t{3}\t{4}\t{5:.2f}'.format(
        n, clid, tally[clid]['seq_num'], tally[clid]['length'], tally[clid]['max_len'], float(tally[clid]['length'])/tally[clid]['seq_num'])


