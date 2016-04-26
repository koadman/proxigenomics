#!/usr/bin/env python
import truthtable as tt
import numpy as np
import sys

table = tt.read_truth(sys.argv[1])
table.print_tally()

