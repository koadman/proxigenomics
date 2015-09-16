#!/usr/bin/env python
import sys
import argparse
from collections import Counter
import pprint as pp

class IntegerSeries:
    '''
    Simple object to issue integer ids. I am sure this exists somewhere in Python
    All series begin at 0
    '''
    def __init__(self, inc=True, step=1):
        '''
        Ascending or descending integers using a given step size
        :param inc: True=ascending, False=descending
        :param step: size of increment
        '''
        self.current = 0
        self.inc = inc
        self.step = step

    @property
    def next(self):
        '''
        Generate the next id and return it.
        :return: next id in series
        '''
        if self.inc:
            self.current += self.step
        else:
            self.current -= self.step
        return self.current


def metisCl_to_mcl(node_id_to_name, cluster_file, output):

    node_id_to_name = {}
    for line in args.node_map[0]:
        cl_assigned = line.rstrip().split()
        if len(cl_assigned) != 2:
            sys.stderr.write('Error: invalid line did not contain 2 fields [{0}]\n'.format(line))
        node_id_to_name[int(cl_assigned[0])] = cl_assigned[1]

    # create a dictionary of scaffold => cluster
    cl_table = Counter()
    iso_id = IntegerSeries(inc=False)
    node_name_to_cl = {}
    for n, line in enumerate(args.cluster_file[0], start=1):
        cl_assigned = [int(ti) for ti in line.rstrip().split()]
        if len(cl_assigned) == 0:
            cl_assigned = [iso_id.next]
            if iso_id.current in cl_table:
                raise RuntimeError('Generated unique id:{0} for nid:{1} already existed'.format(n, iso_id.current))
        node_name_to_cl[node_id_to_name[n]] = cl_assigned
        cl_table.update(cl_assigned)

    print '\nClustering breakdown'
    print 'BEGIN_TABLE'
    print 'cluster\tcount'
    for clid, count in cl_table.items():
        print '{0}\t{1}'.format(clid, count)
    print 'END_TABLE\n'

    # invert the dictionary, to cluster => scaffold
    cl_to_node = dict.fromkeys(cl_table)
    for k, v in node_name_to_cl.iteritems():
        for clid in v:
            cl = cl_to_node.get(clid)
            if cl is None:
                cl = list()
                cl_to_node[clid] = cl
            cl.append(k)

    # sort clusters by descending size
    cl2scf = sorted(cl_to_node.values(), key=len, reverse=True)

    # write to output stream
    for i in range(0, len(cl2scf)):
        args.output.write(' '.join(cl2scf[i]) + '\n')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Convert Metis-style clustering file to MCL format')
    parser.add_argument('node_map', type=argparse.FileType('r'), nargs=1,
                        help='node-map file, implicit ids to names')
    parser.add_argument('cluster_file', type=argparse.FileType('r'), nargs=1,
                        help='Clustering result file')
    parser.add_argument('output', type=argparse.FileType('w'), nargs='?', default=sys.stdout,
                        help='Output file')
    args = parser.parse_args()

    metisCl_to_mcl(*vars(args))