#!/usr/bin/env python
import sys, traceback

import pysam
import numpy
from scipy.stats import histogram
import matplotlib.pyplot as plt


if len(sys.argv) != 2:
    print 'Usage: [bam file]'
    sys.exit(1)

samfile = pysam.AlignmentFile(sys.argv[1], 'rb')

info = []
counts = {'sec': 0, 'unmap': 0, 'rev': 0, 'dup': 0, 'qcf': 0}

for read in samfile.fetch():
    try:
        if read.cigarstring is not None:
            match_len = sum([cig_i[1] for cig_i in read.cigartuples if cig_i[0] == 0])
            read_len = sum([cig_i[1] for cig_i in read.cigartuples])
        else:
            match_len, read_len = (None, None)

        counts['sec'] += int(read.is_secondary)
        counts['unmap'] += int(read.is_unmapped)
        counts['rev'] += int(read.is_reverse)
        counts['dup'] += int(read.is_duplicate)
        counts['qcf'] += int(read.is_qcfail)

        info.append((read.mapping_quality, match_len, read_len))

    except TypeError as e:
        print 'exception {0}'.format(e)
        traceback.print_exc(file=sys.stdout)
        sys.exit(1)

info = numpy.array(info, dtype=numpy.float)

print 'Counts'
print '\taligned reads= {0}\n\tsecondary= {1:.0f}\n\tunmapped= {2:.0f}\n\tdup= {2:.0f}\n\tqcfail= {2:.0f}'.format(
    len(info), counts['sec'], counts['unmap'], counts['dup'], counts['qcf'])

dat = info[:, 0]
dat = dat[~numpy.isnan(dat)]
print 'MapQ stats\n\tmean= {0:.2f} sd= {1:.2f}'.format(numpy.mean(dat), numpy.std(dat))
print '\t50th percentile mapQ= {0}'.format(numpy.percentile(dat, q=50))

#plt.figure(1)
#plt.subplot(211)
#plt.title('Mapping quality distribution')
#plt.hist(dat, bins=20, range=(0, 60))

dat = info[:, 1:3]
dat = dat[~numpy.isnan(dat).any(1)]
rellen = dat[:, 0] / dat[:, 1]

print 'Matching length\n\tmean= {0:.2f} sd= {1:.2f}'.format(numpy.mean(dat[:, 0]), numpy.std(dat[:, 0]))
print 'Fully aligning reads\n\tcount= {0}\n\t50th percentile= {1}'.format(numpy.sum(rellen >= 1.0),
                                                                    numpy.percentile(dat[:, 0], q=50))
#plt.subplot(212)
#plt.title('Distribution of relative mapped read length')
#plt.hist(rellen, bins=20, range=(0, 1))
#plt.show()
