from hic import *
import pysam
import networkx as nx
import sys
import math

if len(sys.argv) != 4:
    print 'Usage: [R1 bam] [R2 bam] [graphml out]'
    sys.exit(1)

# Open the BAM files for each read set.
samR1 = pysam.Samfile(sys.argv[1], 'rb')
samR2 = pysam.Samfile(sys.argv[2], 'rb')

# Build an unconnected graph with scaffolds as vertexes
print 'Building nodes from scaffold names'
if set(samR1.references) != set(samR2.references):
    raise Exception('Error: reference sequences in read-sets are not the same')

g = nx.Graph()
for idx in range(len(samR1.references)):
    g.add_node(samR1.references[idx], length=samR1.lengths[idx])
print "Created {0} nodes".format(g.number_of_nodes())

fragmentRegistry = RegisteredObjectFactory(Fragment)

print 'Scanning R1 for placements'
#fragments = {}
for ar in samR1.fetch():

    # only consider primary alignments
    if ar.is_secondary or ar.is_unmapped:
        continue

    # Note, using a set here eliminates self-loops
    scfName = samR1.getrname(ar.tid)
    frg = fragmentRegistry.requestObject(name=ar.qname)
    if frg.read1:
        raise Exception('Fragment {0} has degenerate placements of R1'.format(frg))
    frg.read1 = ReadPlacement(scfName, ar.qstart)
    frg.read1.mapq = ar.mapq


#print '{0} potential linkages'.format(len(fragments))
print '{0} potential linkages'.format(len(fragmentRegistry))

print 'Scanning R2 for placements'
for ar in samR2.fetch():

    # only consider primary alignments
    if ar.is_secondary or ar.is_unmapped:
        continue

    # record, possibly the second end of fragment.
    scfName = samR2.getrname(ar.tid)
    frg = fragmentRegistry.requestObject(name=ar.qname)
    if frg.read2:
        raise Exception('Fragment {0} has degenerate placements of R2'.format(frg))
    frg.read2 = ReadPlacement(scfName, ar.qstart)
    frg.read2.mapq = ar.mapq

samR1.close()
samR2.close()


def phredq_to_prob(qual):
    return math.pow(10.0, -0.1*qual)


def prob_to_phredq(prob, maxq=None):
    qual = -10.0 * math.log10(prob)
    if maxq and qual > maxq:
        return maxq
    return qual

print 'Adding edges to graph (ignoring self-loops)'
edges = 0
for frg in fragmentRegistry.elements():
    if frg.isInterContig():
        u = frg.read1.contig
        v = frg.read2.contig

        # use conditional prob of r1,r2 map qualities as an edge quality
        # where we cap the value at 60
        p1 = phredq_to_prob(frg.read1.mapq)
        p2 = phredq_to_prob(frg.read2.mapq)
        mapq = prob_to_phredq(p1 * p2, 60)

        if g.has_edge(u, v):
            # weight is a count of coincident edges
            g[u][v]['weight'] += 1
            # sum the mapq score for each coincident
            g[u][v]['mapq'] += mapq
        else:
            g.add_edge(u, v, weight=1, mapq=mapq)
            edges += 1

# traverse edges, convert summed mapq per edge to mean mapq
# over observed coincident edges.
for u, v, dat in g.edges_iter(data=True):
    avg_mapq = int(dat['mapq'] / dat['weight'])
    g[u][v]['mapq'] = avg_mapq

print 'Created {0} inter-scaffold edges'.format(edges)

print 'Writing graphml'
nx.write_graphml(g, sys.argv[3])
nx.write_weighted_edgelist(g, 'edges.csv')
