#!/usr/bin/env python
import truthtable as tt
import numpy as np
import sys

table = tt.read_truth(sys.argv[1])
num_class = []
for k, v in table.asgn_dict.iteritems():
    num_class.append(len(v))

print '#objects\tmean_size\tsd_size'
print '{0}\t{1}\t{2}'.format(len(num_class), np.mean(num_class), np.std(num_class))

