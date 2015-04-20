#!/usr/bin/env python
import numpy as np
import networkx as nx
import argparse


def approx_intra_sim(v, g):
    v_star = g[v]
    deg_vstar = len(v_star)
    if deg_vstar < 1:
        raise RuntimeError('Node degree is than one for node {0}'.format(v))
    # denominator in paper is |v|-1 but doing so means we must exclude all
    # degree 1 nodes (not just isolated |v|=0 nodes)
    return np.sum([dat['weight'] for u1, u2, dat in g.edges(v, data=True) if u2 is not v]) / float(deg_vstar)


def relative_density(v, g):
    v_star = g[v]
    deg_v = len(v_star)
    density_v = len([u for u in v_star if len(g[u]) < deg_v])
    return density_v / float(deg_v)


def relative_compactness(v, g):
    v_star = g[v]
    deg_v = len(v_star)
    sim_v = approx_intra_sim(v, g)
    compact_v = len([u for u in v_star if u is not v and sim_v >= approx_intra_sim(u, g)])
    return compact_v / float(deg_v)


def relevance(v, g):
    return 0.5 * (relative_density(v, g) + relative_compactness(v, g))


class Cover:
    def __init__(self, g):
        self.g = g
        self.hubs = nx.Graph()

    def update(self, vlist):
        for v in vlist:
            self.add(v)

    def add(self, v):
        """
        For the given vertex, build its star graph and add this
        to the list of covering stars for graph.
        :param v: the hub vertex
        """

        # create new node or update existing one
        if v not in self.hubs:
            self.hubs.add_node(v, attr_dict={'analyzed': False, 'hub': True, 'owned_by': set()})
        self.hubs.node[v]['owned_by'].add(v)
        self.hubs.node[v]['hub'] = True

        v_adj = g[v].keys()
        for vi in v_adj:
            # for each satellite, create new node or update existing.
            # we might be covering the same node with more than one star
            if vi not in self.hubs:
                self.hubs.add_node(vi, attr_dict={'analyzed': False, 'hub': False, 'owned_by': set()})
            self.hubs.node[vi]['owned_by'].add(v)
            self.hubs.add_edge(v, vi)

    def shared(self, v):
        """
        Return the list of adjacent vertices (satellites) which are also shared with other hubs
        in cover.
        :param v: the hub vertex
        :return: the list of shared satellites
        """
        v_adj = self.g[v]
        # vertices with more than one owner
        return [u for u in v_adj if len(self.hubs.node[u]['owned_by']) > 1]

    def unshared(self, v):
        """
        Return the list of adjacent vertices (satellites) which are not covered by any other
        hub in cover.
        :param v: the hub vertex
        :return: the list of unshared satellites
        """
        v_adj = self.g[v]
        # vertices with only 1 owner
        return [u for u in v_adj if len(self.hubs.node[u]['owned_by']) <= 1]

    def is_useful(self, v):
        """
        A useful hub is one which has more unshared than shared satellites in the
        cover of G.
        :param v: the hub vertex to test.
        :return: True if this hub is a useful star.
        """
        ns = len(self.shared(v))
        nu = len(self.unshared(v))
        return ns <= nu

    def highest_containing(self, u):
        """
        For a given vertex u, find which owning hub has the highest degree. If there is
        tie between owning hubs, pick one at random.
        :param u: the vertex to test.
        :return: the hub vertex with high degree
        """
        owners = self.hubs.node[u]['owned_by']

        # owners degrees, not include self
        owner_deg = dict((v, self.g.degree(v)) for v in owners if v != u)

        if len(owner_deg) == 0:
            raise RuntimeWarning('{0} is an isolated vertex containing only itself'.format(u))
        elif len(owner_deg) == 1:
            # only one non-self owner
            return owner_deg.popitem()
        else:
            # sort owners by descending degree
            deg_desc = sorted(owner_deg, key=owner_deg.get, reverse=True)
            top = [owner_deg[deg_desc[0]], owner_deg[deg_desc[1]]]
            if top[0] == top[1]:
                # in event of tie, pick winner at random
                idx = np.random.choice(deg_desc[:2])
                return idx, owner_deg[idx]
            else:
                return deg_desc[0], top[0]

    def demote(self, v):
        """
        Demote a node from being a hub. This means all subordinate nodes
        which recognised this node as a hub, are also updated.
        :param v: the target hub to remove
        """
        # for adjacent satellites vi of hub v
        v_adj = set(self.hubs[v])
        for vi in v_adj:
            vi_dat = self.hubs.node[vi]
            if v in vi_dat['owned_by']:
                # TODO this may not be necessary, the logic should result in statement always being true
                vi_dat['owned_by'].remove(v)
            self.hubs.remove_edge(v, vi)

    def migrate(self, vlist, from_v, to_v):
        """
        Update the edges and attributes of a list of nodes to
        reflect the change of ownership.
        :param vlist: the list of target vertices to move
        :param from_v: the origin hub vertex
        :param to_v: the destination hub vertex
        :return:
        """
        for vi in vlist:
            # avoid moving self
            if vi != from_v:
                self.hubs.remove_edge(from_v, vi)
                self.hubs.add_edge(to_v, vi)
                self.hubs.node[vi]['owned_by'].remove(from_v)
                self.hubs.node[vi]['owned_by'].add(to_v)

    def is_hub(self, v):
        return self.hubs.node[v]['hub']

    def is_analyzed(self, v):
        return self.hubs.node[v]['analyzed']

    def mark_analyzed(self, v):
        self.hubs.node[v]['analyzed'] = True

    def exists(self, v):
        """
        The vertex v exists in cover graph.
        :param v: the target vertex.
        :return: True if the vertex exists in the cover graph.
        """
        return v in self.hubs

    def contains(self, vset):
        """
        The set of vertices exists in the cover graph.
        :param vset: the target set of vertices.
        :return: True if the target set is contained by the cover graph.
        """
        return vset < set(self.hubs.nodes())

    def write(self, path):
        gout = self.hubs.copy()
        for v in gout:
            # TODO
            # networkx cannot write collections as attributes
            # so we're converting to a string for now.
            gout.node[v]['owned_by'] = ' '.join(gout.node[v]['owned_by'])
        nx.write_graphml(gout, path)

    def __str__(self):
        return str(self.hubs)

    def __repr__(self):
        return repr(self.hubs)


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Cluster similarity graph using OClustR')
    parser.add_argument('--debug', action='store_true', default=False, help='Enable debugging')
    parser.add_argument('--no-isolates', action='store_true', default=False, help='Ignore isolates in output')
    parser.add_argument('graph', nargs=1, help='Input graphml file')
    parser.add_argument('output', nargs=1, help='Output base file name')
    args = parser.parse_args()

    # read the graph
    g = nx.read_graphml(args.graph[0])
    if args.debug:
        print 'Graph vertices {0} edges {1}'.format(g.order(), g.size())

    # calculate relevance for all vertices
    rel_dict = dict((v, relevance(v, g)) for v in g.nodes() if len(g[v]) >= 1)

    # sort candidates by decreasing relevance
    desc_relevance = sorted(dict((k, v) for k, v in rel_dict.iteritems() if v > 0.0), key=rel_dict.get, reverse=True)

    if args.debug:
        print 'Vertices sorted by decreasing relevance'
        for n, li in enumerate(desc_relevance, start=1):
            print '{0} {1} {2:.3f}'.format(n, li, rel_dict[li])

    # add all isolated nodes
    cover = Cover(g)
    cover.update(set(v for v in g.nodes_iter() if g.degree(v) == 0))
    print 'There were {0} isolated vertices'.format(len(cover.hubs))

    # traverse the set of vertices in descending order of relevance
    # looking for stars which cover new territory
    n_hub = 0
    for v in desc_relevance:
        if not cover.exists(v):
            # new hub of a star
            cover.add(v)
            n_hub += 1
        elif not cover.contains(set(g[v].keys())):
            # star had at least one new covered vertex
            cover.add(v)
            n_hub += 1

    print 'Finished step 1: {0} stars and {1} vertices'.format(n_hub, len(cover.hubs))

    # determine degree of each star and sort by descending degree
    deg_dict = dict((v, g.degree(v)) for v in cover.hubs)
    desc_degree = sorted(deg_dict, key=deg_dict.get, reverse=True)

    if args.debug:
        for n, v in enumerate(desc_degree, start=1):
            if v in rel_dict:
                print '{0} {1} {2:.3f} {3}'.format(n, v, rel_dict[v], deg_dict[v])
            else:
                print '{0} {1} NA {2}'.format(n, v, deg_dict[v])

    # traverse the set of hubs by descending degree, if not useful, then
    # rehome unshared satellites with highest degree overlapping hub.
    demoted = []
    migrated = []
    # for each hub by descending degree
    for v in desc_degree:
        # get the members
        v_adj = g[v]
        for u in v_adj:
            if cover.is_hub(u) and not cover.is_analyzed(u):
                if not cover.is_useful(u):
                    try:
                        v_high, v_deg = cover.highest_containing(u)
                        unshared = cover.unshared(u)
                        if len(unshared) > 0:
                            if args.debug:
                                print 'Moving unshared vertices {0} from {1} to {2}'.format(unshared, u, v_high)
                            cover.migrate(unshared, u, v_high)
                            migrated.extend(unshared)
                        cover.demote(u)
                        demoted.append(u)
                    except RuntimeWarning as wrn:
                        if args.debug:
                            print 'No action taken, {0}'.format(wrn.message)
                cover.mark_analyzed(u)

    print 'Demoted {0} hub vertices, migrated {1} satellites'.format(len(demoted), len(migrated))
    if args.debug:
        print 'Demoted hubs were: {0}'.format(demoted)
        print 'Migrated vertices were: {0}'.format(migrated)

    # write cover graph to file
    cover.write(args.output[0] + '.graphml')

    # write MCL format cluster solution
    with open(args.output[0] + '.mcl', 'w') as out:
        min_deg = 1 if args.no_isolates else 0
        # build a list of hubs and their degree, sort be descending degree for output
        hubs = cover.hubs
        hub_dict = dict((v, hubs.degree(v)) for v in hubs if hubs.node[v]['hub'] and hubs.degree(v) >= min_deg)
        desc_deg = sorted(hub_dict, key=hub_dict.get, reverse=True)
        for vi in desc_deg:
            out.write('{0} {1}\n'.format(vi, ' '.join(cover.hubs[vi])))
        print 'Final solution: {0} stars and {1} satellites'.format(len(hub_dict), sum(hub_dict.values()))
