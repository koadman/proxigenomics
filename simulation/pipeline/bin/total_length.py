#!/usr/bin/env python
import dendropy as dpimport sysdef edge_sum(t):    return sum([ei.length for ei in t.edges() if ei.length])if len(sys.argv) != 2:    print 'Usage: [newick tree]'    sys.exit(1)t = dp.Tree.get(path=sys.argv[1], schema='newick')print edge_sum(t)

