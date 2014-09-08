import vcf

#
# Simple method to count the number of records in a
# VCF file. This is useful for providing progress on
# very large runs.
#
def count_vcf_records(fname):
    return len([rec for rec in vcf.Reader(filename=fname)])


#
# Represents the placement of a read within an
# assembly or mapping.
#
class ReadPlacement:
    def __init__(self, contig, position):
        self.contig = contig
        self.position = position
        self.snpSet = set()
        self.snpInstances = {}

    def __hash__(self):
        return hash((self.contig, self.position))

    def __eq__(self, other):
        if type(other) is not type(self):
            return NotImplemented
        return other.contig == self.contig and other.position == self.position

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return 'SCF_NAME:{0} SCF_START:{1} SNP_INST:{2}'.format(self.contig, self.position, str(self.snpInstances))

    def __repr__(self):
        return self.contig + '-' + str(self.position)

    def addSnp(self, snp):
        if snp in self.snpSet:
            raise Exception('SNP:{0} already exists for READ:{1}'.format(str(snp), str(self)))
        self.snpSet.add(snp)

    def addSnpInstance(self, snp, base):
        if snp in self.snpInstances:
            raise Exception('SNP:{0} already exists as in instance base:{1} for READ:{2}'.format(
                str(snp), self.snpInstances[snp], str(self)))
        self.snpInstances[snp] = base

    def snpCount(self):
        return len(self.snpSet)

#
# Represents a fragment generated from HiC proximity ligation.
# Effectively, this is the same as any other DNA insert.
#
class Fragment:

    def __init__(self, name):
        self.name = name
        self.read1 = None
        self.read2 = None

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        if type(other) is not type(self):
            return NotImplemented
        return other.name == self.name

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return '{0} R1:{1} R2:{2}'.format(self.name, self.read1, self.read2)

    def __repr__(self):
        return self.name

    def isPaired(self):
        return self.read1 is not None and self.read2 is not None

    def isInterContig(self):
        return self.isPaired() and self.read1.contig != self.read2.contig

#
# A SNP site, defined by contig and position
#
class SNP:

    def __init__(self, vcfRecord):
        self.vcfRecord = vcfRecord
        self.reference = str(vcfRecord.REF).upper()
        if len(vcfRecord.ALT) != 1:
            raise Exception('Multi-allelic variants are not supported: {0}'.format(vcfRecord))
        self.variant = str(vcfRecord.ALT[0]).upper()
        self.alleles = {}

    @property
    def contig(self):
        return self.vcfRecord.CHROM

    @property
    def position(self):
        return self.vcfRecord.POS

    @property
    def quality(self):
        return float(self.vcfRecord.QUAL)

    @property
    def ref_id(self):
        return '{0}:{1}'.format(self, self.reference)

    @property
    def var_id(self):
        return '{0}:{1}'.format(self, self.variant)

    def get_split_ids(self):
        return [self.ref_id, self.var_id]

    def __eq__(self, other):
        if type(other) is not type(self):
            return NotImplemented
        return other.contig == self.contig and other.position == self.position

    def __ne__(self,other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.contig, self.position))

    def __str__(self):
        return '{0}:{1}'.format(self.contig, self.position)

    def __repr__(self):
        return self.contig + '-' + str(self.position)

    # check that the candidate allele is either the reference or identified
    # variant. This deals with noise in the pileup, where a small rate of aligned
    # reads possess neither the reference or identified variant.
    def isUndefinedAllele(self, allele):
        return not (allele == self.reference or allele == self.variant)

    def addAllele(self, candidate, fragment):
        candidate = candidate.upper()
        if candidate not in self.alleles:
            self.alleles[candidate] = []
        self.alleles[candidate].append(fragment)

    @property
    def weights(self):
        w = {}
        for k, v in self.alleles.iteritems():
            w[k] = len(v)
        return w

    @property
    def ratio(self):
        w = self.weights
        cv = w.get(self.variant, 0.0)
        if cv == 0.0:
            return 0.0
        cr = w.get(self.reference, 0.0)
        return float(cv)/float(cr+cv)

    @property
    def depth(self):
        return sum(self.weights.values())

#
# This could probably be made into a builder.
# - we want to create objects if new or return preexisting instances
#   if the object identity is equal.
#
# With a type and attr, could call type(classname, object, attr*)
#
class RegisteredObjectFactory:

    def __init__(self, clazz):
        self.clazz = clazz
        self.registry = {}

    def __len__(self):
        return len(self.registry)

    #
    # Create a new object instance and register it.
    #
    def requestObject(self, **kwargs):
        obj = self.clazz(**kwargs)
        if obj in self.registry:
            obj = self.registry[obj]
        else:
            self.registry[obj] = obj
        return obj

    #
    # Get the object or raise exception
    #
    def getObject(self, **kwargs):
        obj = self.clazz(**kwargs)
        return self.registry[obj]

    def elements(self):
        return self.registry.values()
