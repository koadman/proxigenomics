from hic import *
import networkx as nx
import datetime as dt
import os.path
import argparse
import logging
import pysam
import vcf
import sys

# log config
logging.basicConfig(filename='snpNetwork.log',level=logging.DEBUG)


#
# Read base and quality
#
def read_base_and_quality(alignedRead, pos):
    return alignedRead.seq[pos], ord(alignedRead.qual[pos]) - 33


#
# Test for file existence or raise an error
#
def file_exists(fnames):
    for fn in fnames:
        if not os.path.exists(fn):
            raise IOError('Error: \"{0}\" does not exist'.format(fn))

#
# User interface
#
parser = argparse.ArgumentParser(description='Build snp graph from HiC sequencing data')
parser.add_argument('--split', dest='split_node', help='Binary variant states are treated as independent nodes',
                    action='store_true', default=False)
parser.add_argument('vcf_file', help='VCF file of predicted variant sites', metavar='VCF_FILE')
parser.add_argument('r1_file', help='BAM file for R1 reads', metavar='R1_BAM')
parser.add_argument('r2_file', help='BAM file for R2 reads', metavar='R2_BAM')
parser.add_argument('output', help='Output graph')
args = parser.parse_args()


try:
    # check input files exist
    file_exists([args.vcf_file, args.r1_file, args.r2_file])

    vcfFile = vcf.Reader(filename=args.vcf_file)
    samR1 = pysam.Samfile(args.r1_file, 'rb')
    samR2 = pysam.Samfile(args.r2_file, 'rb')

    varCount = 0
    firstTime = dt.datetime.now()
    lastTime = firstTime

    # Registry objects for tracking fragments and snps
    frgRegistry = RegisteredObjectFactory(Fragment)
    snpRegistry = RegisteredObjectFactory(SNP)

    for variant in vcfFile:

        # Skip any variant that isn't a SNP
        if not variant.is_snp:
            continue

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
                base, bq = read_base_and_quality(aln, pr.qpos)
                mq = aln.mapq

                if bq < 30:
                    continue

#                if aln.qname == 'M00958:3:000000000-A44NP:1:2110:22807:10490':
#                    print 'r1',snpPos, base, bq, mq
#                    print aln

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

                base, bq = read_base_and_quality(aln, pr.qpos)
                mq = aln.mapq

                if bq < 30:
                    continue

#                if aln.qname == 'M00958:3:000000000-A44NP:1:2110:22807:10490':
#                    print 'r2',snpPos, base, bq, mq
#                    print aln

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

    if args.split_node:

        for frg in frgRegistry.elements():

            if not frg.isPaired():
                continue

            print frg, '\n\n'

            for snpR1 in frg.read1.snpSet:

                # create ref node
                # create var node
                print snpR1, '\n\n'
                print snpR1.reference, snpR1.alleles[snpR1.reference], "\n\n"
                print snpR1.variant, snpR1.alleles[snpR1.variant], "\n\n"

                print frg.name in snpR1.alleles[snpR1.reference]
                print frg.name in snpR1.alleles[snpR1.variant]
                sys.exit(1)

                for snpR2 in frg.read2.snpSet:

                    # skip self loops
                    if snpR1 == snpR2:
                        continue


                    for id1 in snpR1.get_split_ids():
                        for id2 in snpR2.get_split_ids():
                            if g.has_edge(id1, id2):
                                g[id1][id2]['weight'] += 1
                            else:
                                g.add_edge(id1, id2, weight=1)

    else:

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

    nx.write_graphml(g, args.output)

except Exception as e:
    print e
    sys.exit(1)
