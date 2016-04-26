#!/usr/bin/env python
import networkx as nx
import community as com
import numpy as np
import truthtable as tt

def decompose_graph(g):
    """
    Using the Louvain algorithm for community detection, as
    implemented in the community module, determine the partitioning
    which maximises the modularity. For each individual partition
    create the sub-graph of g

    :param g: the graph to decompose
    :return: the set of sub-graphs which form the best partitioning of g
    """
    decomposed = []
    p = com.best_partition(g)
    partition_labels = np.unique(p.values())

    # for each partition, create the sub-graph
    for pi in partition_labels:
        # start with a complete copy of the graph
        gi = g.copy()
        # build the list of nodes not in this partition and remove them
        to_remove = [n for n in g.nodes_iter() if p[n] != pi]
        gi.remove_nodes_from(to_remove)
        decomposed.append(gi)

    return decomposed


def cluster(g, no_iso, soft=True):

    # remove isolated nodes
    if no_iso:
        singletons = [n for n, degree in g.degree().items() if degree < 1]
        print 'Removed {0} isolated nodes from graph'.format(len(singletons))
        g.remove_nodes_from(singletons)
        print_info(g)

    subgraphs = list(nx.connected_component_subgraphs(g))
    print 'The graph comprises {0} connected components'.format(len(subgraphs))
    print 'Components with more than one node:'
    for n, sg in enumerate(subgraphs, start=1):
        if sg.order() > 1:
            print '\tcomponent {0}: {1} nodes {2} edges'.format(n, sg.order(), sg.size())

    # determine the best partitioning of g
    partitions = com.best_partition(g)

    # build a dictionary of classification from the partitioning result
    # this is effectively a hard clustering answer to the problem
    com_ids = set(partitions.values())
    print 'There were {0} communities'.format(len(com_ids))
    print 'Communities with more than one node:'

    communities = {}
    for ci in com_ids:
        communities[ci] = dict((n, 1.0) for n, cj in partitions.iteritems() if cj == ci)
        if len(communities[ci]) > 1:
            print '\tcommunity {0}: {1} nodes'.format(ci, len(communities[ci]))

    if soft:
        # iterate over the nodes in the graph, if an edge connects nodes in disjoint
        # partitions, then make both nodes members of both partitions.
        for n1 in g.nodes_iter():
            for n2 in nx.all_neighbors(g, n1):
                if partitions[n1] != partitions[n2]:
                    # we add them at half weight just for discrimination
                    communities[partitions[n1]][n2] = 0.5
                    communities[partitions[n2]][n1] = 0.5

    return communities


def print_info(g):
    """
    Print order and size of graph g
    :param g: graph to print info
    """
    print 'Graph composed of {0} nodes and {1} edges'.format(g.order(), g.size())


def write_mcl(communities, path):
    """
    Write communities dictionary to output file.
    :param communities: the dictionary of communities
    :param path: the output file path
    """
    with open(path, 'w') as hout:
        keys = sorted(communities.keys())
        for k in keys:
            line = ' '.join(sorted(communities[k].keys()))
            hout.write(line.strip())
            hout.write('\n')


def plot_graph(part, g):
    import matplotlib.pyplot as plt
    from palettable.colorbrewer.qualitative import Paired_12, Set3_12

    # use 24 colours for plotting partitions.
    color_list = Paired_12.mpl_colors + Set3_12.mpl_colors
    ncol = len(color_list)

    layouts = [nx.graphviz_layout, nx.spring_layout]

    for i in xrange(2):
        pos = layouts[i](g)
        plt.subplot(121 + i)
        for n, com in enumerate(set(part.values())):
            list_nodes = [nodes for nodes in part.keys() if part[nodes] == com]
            nx.draw_networkx_nodes(g, pos, list_nodes, node_size=30, node_color=color_list[n % ncol], alpha=1.0)
        nx.draw_networkx_edges(g, pos, alpha=0.5)
    plt.show()


def write_output(communities, filename, ofmt='mcl'):
    if ofmt == 'mcl':
        write_mcl(communities, filename)
    elif ofmt == 'graphml':
        # convert communities to a graph
        cg = nx.DiGraph()
        for k, v in communities.iteritems():
            cg.add_node(k)
            for vi in v:
                cg.add_edge(k, vi)
        nx.write_graphml(cg, filename)
    else:
        raise RuntimeError('Unsupported format type: {0}'.format(ofmt))

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Decompose a graph into its communities')
    parser.add_argument('--no-isolates', action='store_true', default=False, help='Remove isolated nodes')
    parser.add_argument('--otype', choices=['hard', 'soft', 'induced'], default='soft', help='Output type')
    parser.add_argument('--ofmt', choices=['mcl', 'graphml'], default='mcl', help='Specify output format [mcl]')
    parser.add_argument('-b', '--threshold', type=float, help='Edge weight threshold to apply before clustering')
    parser.add_argument('input', nargs=1, help='Input graph (graphml format)')
    parser.add_argument('output', nargs=1, help='Output file')
    args = parser.parse_args()

    if args.otype == 'induced':
        raise RuntimeError('induced option no longer supported')

    g = nx.read_graphml(args.input[0])
    print 'Initial statistics'
    print_info(g)

    if args.threshold:
        thres_g = g.copy()
        # thresholding, current normalises first because we've got raw weights
        for u,v,d in g.edges_iter(data=True):
            thres_g[u][v]['weight'] = float(d['weight']) / (float(g.node[u]['length'] * g.node[v]['length']) / 65536.0)
            if thres_g[u][v]['weight'] < args.threshold:
                thres_g.remove_edge(u, v)
        with open('test.mcl', 'w') as out_h:
            for u,v,d in thres_g.edges_iter(data=True):
                out_h.write('{0} {1} {2}\n'.format(u, v, d['weight']))
        print 'Thresholding reduced graph edge count from {0} to {1}'.format(g.size(), thres_g.size())
        g = thres_g

    communities = cluster(g, args.no_isolates, args.otype == 'soft')

    write_output(communities, args.output[0], args.ofmt)

    # Write out communities as a hard clustering solution
    #if args.otype == 'hard':
    #        write_mcl(communities, args.output[0])

    #elif args.otype == 'induced':
        # Find nodes which connect to multiple communities.
        # Add their membership to each and use this as a soft
        # clustering solution.
        # induced = com.induced_graph(partitions, g)
        # sub_ind = nx.connected_component_subgraphs(induced)
        # groups = {}
        # n = 1
        # for sg in sub_ind:
        #
        #     # skip communities of one node
        #     if sg.size() <= 2:
        #         continue
        #
        #     cut, parts = nx.stoer_wagner(sg)
        #     for p in parts:
        #         if n not in groups:
        #             groups[n] = {}
        #         for pi in p:
        #                 dt = dict((nid, 1.0) for nid, cj in partitions.iteritems() if cj == pi)
        #                 groups[n].update(dt)
        #         n += 1
        #
        # write_mcl(groups, args.output[0])

    #elif args.otype == 'soft':
    #    # iterate over the nodes in the graph, if an edge connects nodes in disjoint
    #    # partitions, then make both nodes members of both partitions.
    #    for n1 in g.nodes_iter():
    #        for n2 in nx.all_neighbors(g, n1):
    #            if partitions[n1] != partitions[n2]:
    #                # we add them at half weight just for discrimination
    #                communities[partitions[n1]][n2] = 0.5
    #                communities[partitions[n2]][n1] = 0.5

    #    write_mcl(communities, args.output[0])

    #else:
    #    raise RuntimeError('Unsupported output type')
