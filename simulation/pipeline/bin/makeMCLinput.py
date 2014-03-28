#!/usr/bin/env python
import pandas as pd
import sys
import os

if len(sys.argv) != 4:
	print 'Usage: [min length] [edges csv] [nodes csv]'
	sys.exit(1)

# Minimum sequence length
minLength = int(sys.argv[1])

# Load edge and node tables
edges = pd.read_csv(sys.argv[2],sep=' ')
nodes = pd.read_csv(sys.argv[3],sep=' ')

# Filter out edges that aren't in 'keepers'
keepers = nodes[nodes.LENGTH > minLength]
filtered = edges[(edges.TARGET.isin(keepers.ID) & edges.SOURCE.isin(keepers.ID))]

# Write mcl input files
filtered.to_csv('mclIn.weighted', cols=['SOURCE','TARGET','WEIGHT'], sep=' ', index=False, header=False)
filtered.to_csv('mclIn.unweighted', cols=['SOURCE','TARGET','RAWWEIGHT'], sep=' ', index=False, header=False)
