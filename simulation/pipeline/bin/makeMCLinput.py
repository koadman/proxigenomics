import pandas as pd
import sys
import os

if len(sys.argv) != 3:
	print 'Usage: [min length] [run folder]'
	sys.exit(1)

minLength = int(sys.argv[1])
nodes = pd.read_csv(os.path.join(sys.argv[2],'nodes.csv'),sep=' ')
edges = pd.read_csv(os.path.join(sys.argv[2],'edges.csv'),sep=' ')

# Filter out edges that aren't in 'keepers'
keepers = nodes[nodes.LENGTH > minLength]
filtered = edges[(edges.TARGET.isin(keepers.ID) & edges.SOURCE.isin(keepers.ID))]

# Write mcl input files
filtered.to_csv('mclIn.weighted', cols=['SOURCE','TARGET','WEIGHT'], sep=' ', index=False, header=False)
filtered.to_csv('mclIn.unweighted', cols=['SOURCE','TARGET','RAWWEIGHT'], sep=' ', index=False, header=False)
