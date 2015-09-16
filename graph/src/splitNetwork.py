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
logging.basicConfig(filename='split-network.log', level=logging.DEBUG)


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
# Conditionally add new nodes to the graph, where kwargs contains additional
# key/value attributes for the node.
#
def add_node(g, id, **kwargs):
    if not g.has_node(id):
        g.add_node(id, kwargs)


#
# User interface
#
parser = argparse.ArgumentParser(description='Build snp graph from HiC sequencing data')
parser.add_argument('-b', '--base_quality', help='Minimum base quality', type=int, default=0)
parser.add_argument('-m', '--map_quality', help='Minimum mapping quality', type=int, default=0)
parser.add_argument('-v', '--variant_quality', help='Minimum mapping quality', type=int, default=0)
parser.add_argument('vcf_file', help='VCF file of predicted variant sites', metavar='VCF_FILE')
parser.add_argument('r1_file', help='BAM file for R1 reads', metavar='R1_BAM')
parser.add_argument('r2_file', help='BAM file for R2 reads', metavar='R2_BAM')
parser.add_argument('output', help='Output graph')
args = parser.parse_args()

try:
    # check input files exist
    file_exists([args.vcf_file, args.r1_file, args.r2_file])

    # get the total number of VCF records for progress reporting
    n_vcf = count_vcf_records(args.vcf_file)

    vcfFile = vcf.Reader(filename=args.vcf_file)
    samR1 = pysam.Samfile(args.r1_file, 'rb')
    samR2 = pysam.Samfile(args.r2_file, 'rb')

    # Some user-feedback progress reporting variables
    varCount = 0
    firstTime = dt.datetime.now()
    lastTime = firstTime

    # Registry objects for tracking fragments and snps
    frgRegistry = RegisteredObjectFactory(Fragment)
    snpRegistry = RegisteredObjectFactory(SNP)


    #
    # Here, we iterate over each variant site and find the pileup column
    # within R1 and then R2 bam files. For each site, we track Fragments and
    # there associated reads, as well as the SNP itself.
    #
    # These build up a registry of fragments and snps, accessed through object
    # identity.
    #
    for variant in vcfFile:

        # Skip any variant that isn't a SNP
        if not variant.is_snp:
            continue

        # impose minimum quality on variants. Quality will depend on tool
        # which predicted variant site.
        if variant.QUAL < args.variant_quality:
            logging.info('%s was below quality threshold: vq=%d', variant, variant.QUAL)
            continue

        # register/get a SNP
        try:
            snp = snpRegistry.requestObject(vcfRecord=variant)
        except Exception as ex:
            logging.info('Skipping: %s', ex)
            continue

        # Iterate over R1 mapping, looking at the pileup at SNP position
        for n, col in enumerate(samR1.pileup(reference=snp.contig, start=snp.position-1, end=snp.position, truncate=True)):

            # probably unnecessary, added while trying to understand under-documented API
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

                # impose minimum quality threshold on basecall and alingment
                if bq < args.base_quality or mq < args.map_quality:
                    logging.info('%s was below quality thresholds: bq=%d, mq=%s', base, bq, mq)
                    continue

                # skip undefined variants
                if snp.isUndefinedAllele(base):
                    logging.info('%s was not ref %s nor alt %s', base, snp.reference, snp.variant)
                    continue

                # Obtain fragment from registry
                frg = frgRegistry.requestObject(name=aln.qname)

                rpl = ReadPlacement(snp.contig, aln.pos)
                if not frg.read1:
                    frg.read1 = rpl
                elif rpl != frg.read1:
                    logging.warn('Tried to assign different read placement r1 [%s] to fragment [%s]', rpl, frg)

                # register allele
                frg.read1.addSnpInstance(snp, base)

        # Now iterate over R2 mapping, looking at the pileup at SNP position
        for n, col in enumerate(samR2.pileup(reference=snp.contig, start=snp.position-1, end=snp.position, truncate=True)):

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

                # impose minimum quality threshold on basecall and alingment
                if bq < args.base_quality or mq < args.map_quality:
                    logging.info('%s was below quality thresholds: bq=%d, mq=%s', base, bq, mq)
                    continue

                # skip undefined variants
                if snp.isUndefinedAllele(base):
                    logging.info('%s was not ref %s nor alt %s', base, snp.reference, snp.variant)
                    continue

                # Obtain fragment from registry
                frg = frgRegistry.requestObject(name=aln.qname)

                rpl = ReadPlacement(snp.contig, aln.pos)
                if not frg.read2:
                    frg.read2 = rpl
                elif rpl != frg.read2:
                    logging.warn('Tried to assign different read placement r2 [%s] to fragment [%s]', rpl, frg)

                # register allele
                frg.read2.addSnpInstance(snp, base)

        varCount += 1
        if varCount % 100 == 0:
            curTime = dt.datetime.now()
            ta = (curTime - lastTime)
            tb = (curTime - firstTime)
            print "... processed {0}/{1} variants in {2} walltime: {3}".format(varCount, n_vcf, ta, tb)
            lastTime = curTime

    print 'Finished reading data, {0} variants'.format(varCount)

    samR1.close()
    samR2.close()

    print 'Registered {0} fragments'.format(len(frgRegistry))

    #
    # Now we just build the graph.
    #
    # - Nodes are uniquely defined SNP as: contig, pos and variant base
    # - Edges are defined by fragment (R1,R2) placement, where edge weight represents
    # accumulated fragment count.
    #

    # counting various things causing rejection
    rejCount = {'contradictory': 0, 'self-loop': 0, 'unpaired': 0}

    # empty graph
    g = nx.Graph(type='split', version=1)

    # iterate over all fragments in registry
    for frg in frgRegistry.elements():

        if not frg.isPaired():
            rejCount['unpaired'] += 1
            continue

        # iterate over all SNPs contained within a fragment
        # nested loops as this is n-way.
        for snpR1, baseR1 in frg.read1.snpInstances.iteritems():

            # TODO: labels should be generated in class, this is really a form of object identity
            u = '{0}.{1}'.format(snpR1, baseR1)  # origin

            for snpR2, baseR2 in frg.read2.snpInstances.iteritems():


                # TODO: ditto above
                v = '{0}.{1}'.format(snpR2, baseR2)  # destination

                # Skipping self loops, often occurring if read pairs overlap
                if snpR1 == snpR2:
                    # just interesting to capture contradictions info.
                    if baseR1 != baseR2:
                        logging.warn('Contradictory bases for R1/R2 {0}/{1} at an overlapping SNP position {2}'.format(
                            baseR1, baseR2, snpR1))
                        rejCount['contradictory'] += 1
                        continue
                    else:
                        logging.warn('Self-loop for R1/R2 {0}/{1} at an overlapping SNP position {2}'.format(
                            baseR1, baseR2, snpR1))
                        rejCount['self-loop'] += 1
                        continue

                #
                # WEIGHTING!
                #
                # Weights by simple occurrence
                #
                if g.has_edge(u, v):
                    # add one to edge
                    g[u][v]['weight'] += 1
                else:
                    # pre-create nodes so we can add some attributes
                    # TODO: add more or more pertinent attributes.
                    add_node(g, u, quality=snpR1.quality)
                    add_node(g, v, quality=snpR2.quality)
                    # new edge
                    g.add_edge(u, v, weight=1)


    print "Created {0} nodes".format(g.number_of_nodes())
    print 'Rejected snp instances: {0}'.format(rejCount)
    nx.write_graphml(g, args.output)

except Exception as e:
    print e
    sys.exit(1)
