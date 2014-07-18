from hic import *
import pysam
import vcf
import networkx as nx
import sys
import datetime as dt
import logging

# log config
logging.basicConfig(filename='snpNetwork.log',level=logging.DEBUG)

#
# Obtain a fragment from the registry. If not found
# register it by fragment name.
#
# def registerFragment(frgName):
#     frg = frgRegistry.get(frgName)
#     if frg is None:
#         frg = Fragment(frgName)
#         frgRegistry[frgName] = frg
#     return frg

#
# Read base and quality
#
def readBaseAndQuality(alignedRead, pos):
    return alignedRead.seq[pos], ord(alignedRead.qual[pos]) - 33

if len(sys.argv) != 5:
    print 'Usage: [VCF] [R1 bam] [R2 bam] [graph out]'
    sys.exit(1)

vcfFile = vcf.Reader(filename=sys.argv[1])
samR1 = pysam.Samfile(sys.argv[2], 'rb')
samR2 = pysam.Samfile(sys.argv[3], 'rb')

varCount = 0
firstTime = dt.datetime.now()
lastTime = firstTime
frgRegistry = RegisteredObjectFactory(Fragment)
snpRegistry = RegisteredObjectFactory(SNP)

for variant in vcfFile:

    if variant.is_snp == False:
        continue  # Skip any variant that isn't a SNP

    snpPos = variant.POS
    scfName = variant.CHROM

    try:
        snp = snpRegistry.requestObject(vcfRecord=variant)
    except Exception as ex:
        logging.info('Skipping: %s', ex)
        continue

    # Iterate over R1 mapping, looking at the pileup at SNP position
    for n, col in enumerate(samR1.pileup(reference=scfName, start=snpPos-1, end=snpPos, truncate=True)):

        if n > 1:
            raise Exception('Pileup is assumed to produced only one exact column')

        for pr in col.pileups:

            aln = pr.alignment

            # skip secondary alignments and indels
            if aln.is_secondary or pr.indel != 0:
                logging.debug('Skipped secondary or indel: %s', aln)
                continue

            # nucleotide at SNP site for this read.
            base, bq = readBaseAndQuality(aln, pr.qpos)
            mq = aln.mapq
            #print base, bq, mq

            # skip undefined variants
            if snp.isUndefinedAllele(base):
                logging.info('%s was not ref %s nor alt %s', base, snp.reference, snp.variant)
                continue

            # Obtain fragment from registry
            frg = frgRegistry.requestObject(name=aln.qname)

            rpl = ReadPlacement(scfName, aln.pos)
            if not frg.read1:
                frg.read1 = rpl
            elif rpl != frg.read1:
                logging.warn('Tried to assign different read placement r1 [%s] to fragment [%s]', rpl, frg)

            # register allele
            snp.addAllele(base, aln.qname)
            try:
                frg.read1.addSnp(snp)
            except Exception as ex:
                logging.warn('%s', ex)

    # Now iterate over R2 mapping, looking at the pileup at SNP position
    for n, col in enumerate(samR2.pileup(reference=scfName, start=snpPos-1, end=snpPos, truncate=True)):

        if n > 1:
            raise Exception('Pileup is assumed to produced only one exact column')

        for pr in col.pileups:

            aln = pr.alignment

            # skip secondary alignments and indels
            if aln.is_secondary or pr.indel != 0:
                logging.debug('Skipped secondary or indel: %s', aln)
                continue

            base, bq = readBaseAndQuality(aln, pr.qpos)
            mq = aln.mapq
            #print base, bq, mq

            # skip undefined variants
            if snp.isUndefinedAllele(base):
                logging.info('%s was not ref %s nor alt %s', base, snp.reference, snp.variant)
                continue

            # Obtain fragment from registry
            frg = frgRegistry.requestObject(name=aln.qname)

            rpl = ReadPlacement(scfName, aln.pos)
            if not frg.read2:
                frg.read2 = rpl
            elif rpl != frg.read2:
                logging.warn('Tried to assign different read placement r2 [%s] to fragment [%s]', rpl, frg)

            # register allele
            snp.addAllele(base, aln.qname)
            try:
                frg.read2.addSnp(snp)
            except Exception as ex:
                logging.warn('%s', ex)

    varCount += 1
    if varCount % 100 == 0:
        curTime = dt.datetime.now()
        ta = (curTime - lastTime)
        tb = (curTime - firstTime)
        print "...processed {0} variants in {1} walltime: {2}".format(varCount, ta, tb)
        lastTime = curTime

print 'Finished reading data, {0} variants'.format(varCount)

samR1.close()
samR2.close()

print 'Registered {0} fragments'.format(len(frgRegistry))

#
# Now below we should just build out the graph.
# Nodes are uniquely defined SNP (chrom + pos)
# Edges are defined by fragment (R1,R2) placement.
# if a fragment does not contain a SNP at each end, it is not included.
#
g = nx.Graph()

# nodes/vertexes
for snp in snpRegistry.elements():
    g.add_node(snp)
    g.node[snp]['ref'] = snp.reference
    g.node[snp]['var'] = snp.variant
    g.node[snp]['ratio'] = snp.ratio
    g.node[snp]['qual'] = snp.quality
    g.node[snp]['depth'] = snp.depth

# links/edges
for frg in frgRegistry.elements():

    if not frg.isPaired():
        continue

    for snpR1 in frg.read1.snpSet:
        #print '{0} {1}'.format(snpR1, snpR1.weights)
        for snpR2 in frg.read2.snpSet:
            if snpR1 == snpR2:
                # skip self loops
                continue
            #print '\t{0} {1}'.format(snpR2, snpR2.weights)
            if g.has_edge(snpR1, snpR2):
                g[snpR1][snpR2]['weight'] += 1
            else:
                g.add_edge(snpR1, snpR2, weight=1)

print "Created {0} nodes".format(g.number_of_nodes())

nx.write_graphml(g, sys.argv[4])
