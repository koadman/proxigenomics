# NOTE: This code is not correct.
# 
# I have put this to the side as considering strandedness was doing my head in.
#
from Bio import SeqIO
from Bio import Alphabet
from Bio.Restriction import *
from collections import OrderedDict
import numpy
import sys

# restriction enzymes used in reaction
cut4 = Restriction.NlaIII
cut6 = Restriction.ClaI

def findRestrictionSites(enzyme, seq):
    """For supplied enzyme, find all restriction sites in a given sequence
    returns list of sites.
    """
    return enzyme.search(seq,linear=False)

def drawDelta(minLength,maxLength):
    """Draws from a distribution only accepting values between min and max."""
    # as this relates to a circular chromosome, the min and max could be considered one value.
    # could do this modulo length of chromosome.
    
    geometricProbability = 1.0e-4
    
    delta = numpy.random.geometric(p=geometricProbability,size=1)
    while (delta < minLength and delta > maxLength):
        delta = numpy.random.geometric(p=geometricProbability,size=1)
    return int(delta)

def makeRead(seq, fwdRead, length):
    read = None
    
    # edge case - return whole sequence of read length > sequence
    if length >= len(seq):
        read = seq;
        if not fwdRead:
            read = read.reverse_complement(id=True,description=True)
        return read
    
    # return forward or reverse read
    if fwdRead:
        read = seq[:length]
    else:
        read = seq[-length:]
        read = read.reverse_complement(id=True,description=True)
        
    return read

def writeReads(handle, sequence, outputFormat):
    """Append sequence to file in fastq format. Currently no effort is spent
    to stream these records"""
    SeqIO.write(sequence, handle, outputFormat)

class Part:
    def __init__(self, seq, pos1, pos2, fwd, repIdx):
        self.seq = seq
        self.pos1 = pos1
        self.pos2 = pos2
        self.fwd = fwd
        self.repIdx = repIdx
    def __repr__(self):
        return repr((self.seq,self.pos1,self.pos2,self.fwd,self.repIdx))
    
class Cell:
    'Represents a cell in the community'
    def __init__(self, name, abundance):
        self.name = name
        self.abundance = float(abundance)
    def __repr__(self):
        return repr((self.name,self.abundance))
    def __str__(self):
        return self.name

class Replicon:
    'Represents a replicon which holds a reference to its containing cell'
    def __init__(self, name, parentCell, sequence):
        self.name = name
        self.parentCell = parentCell
        self.sequence = sequence
        self.cutSites = {}
        self.cutSites['4cut'] = numpy.array(findRestrictionSites(cut4, sequence.seq))
        self.cutSites['6cut'] = numpy.array(findRestrictionSites(cut6, sequence.seq))
    
    def __repr__(self):
        return repr((self.name, self.parentCell, self.sequence))
    
    def __str__(self):
        return str(self.parentCell) + '.' + self.name
    
    def length(self):
        return len(self.sequence)
    
    def subSeq(self, pos1, pos2, fwd):
        """Create a subsequence from replicon where positions are strand relative.
        Those with fwd=False will be reverse complemented.
        """
        # Likely to have a off-by-one error in this code.
        ss = None
        if fwd:
            ss = self.sequence[pos1:pos2]
            ss.description = str(pos1) + "..." + str(pos2) + ":" + str(fwd)
        else:
            ss = self.sequence[pos2:pos1] 
            ss.description = str(pos2) + "..." + str(pos1) + ":" + str(fwd)
            ss = ss.reverse_complement(id=True,description=True)
        return ss
    
    def randomCutSite(self, cutType):
        'Return a uniformly random cut-site'
        cs = self.cutSites[cutType]
        idx = numpy.random.randint(low=0,high=len(cs))
        return cs[idx]
    
    def nearestCutSiteAbove(self, cutType, pos, fwd):
        cs = self.cutSites[cutType]
        idx = numpy.searchsorted(cs, pos)
        if fwd:
            return cs[idx]
        else:
            if idx == 0: return cs[-1]
            else: return cs[idx-1]
        
    def nearestCutSiteByDistance(self, cutType, pos):
        """Find the nearest restriction cut site for the specified cutter type [4cut, 6cut]
        returns genomic position of nearest cut site"""
        cs = self.cutSites[cutType]
        idx = numpy.searchsorted(cs, pos)
        d1 = None
        d2 = None
        # edge cases when insertion point between be before or after 
        # beginning or end of linear sequence
        if idx == len(cs):
            d1 = (idx-1, pos - cs[-1])
            d2 = (0, cs[0])
        elif idx == 0:
            d1 = (idx-1, self.length() - cs[-1]) 
            d2 = (0, cs[0])
        else:
            d1 = (idx, pos - cs[idx])
            d2 = (idx+1, cs[idx+1] - pos)

        if d2[1] < d1[1]: return cs[d2[0]]
        else: return cs[d1[0]]

class Community:
    """Represents the community, maintaining a registry of cells and replicons.
    
    Both cells and replicons are expected to be uniquely named across the commmunity. This constraint comes
    from the requirement that a multi-fasta file contain unique sequence identifiers.
    
    A table is used for a basic community definition, with the following three column format. Lines beginning
    with # and empty lines are ignored.
    
    [replicon name] [cell name] [abundance]
    """
    def __init__(self, interRepliconProbability, tableFileName, seqFileName):
        self.createdReadCount = 0
        self.pdf = None
        self.cdf = None
        self.totalRawAbundance = 0
        self.repliconRegistry = OrderedDict()
        self.indexToName = None
        self.cellRegistry = OrderedDict()
        self.interRepliconProbability = interRepliconProbability 
        # Read in the sequences
        sequences = SeqIO.to_dict(SeqIO.parse(open(seqFileName), 'fasta', Alphabet.generic_dna))
        hTable = open(tableFileName,'r')
        for line in hTable:
            line = line.rstrip().lstrip()
            if line.startswith('#') or len(line) == 0:
                continue 
            field = line.split()
            if len(field) != 3:
                print 'sequence table has missing fields at [',  line, ']'
                sys.exit(1)
            repliconName = field[0]
            cellName = field[1]
            cellAbundance = field[2]
            parentCell = self.registerCell(cellName,cellAbundance)
            self.registerReplicon(repliconName, parentCell, sequences.get(repliconName))
        hTable.close()
        self.__initProbabilities()

    
    def incCreatedReadCount(self):
        self.createdReadCount += 1
        
    def registerReplicon(self, name, parentCell, sequence):
        'Add a new cell to the community cell registry'
        replicon = self.repliconRegistry.get(name)
        if replicon == None:
            replicon = Replicon(name,parentCell,sequence)
            self.repliconRegistry[name] = replicon
        return replicon

    def registerCell(self, cellName, cellAbundance):
        'Add a new replicon to the community replicon registry'
        parentCell = self.cellRegistry.get(cellName)
        if parentCell == None:
            parentCell = Cell(cellName,cellAbundance)
            self.cellRegistry[cellName] = parentCell
            self.__updateTotalAbundance()  
        return parentCell
    
    # Total abundance of all cells in registry.
    def __updateTotalAbundance(self):
        'Recalculate the total relative abundance specified for the community by referring to the registry'
        ab = 0
        for ca in self.cellRegistry.values():
            ab += ca.abundance
        self.totalRawAbundance = ab
    
    def getReplicon(self,name):
        'Return a replicon instance by name from the registry'
        return self.__getitem__(name)
    
    def getRepliconLength(self,name):
        'Return the length in basepairs of a replicon by name from the registry'
        return self.repliconRegistry[name].length()
    
    def __initProbabilities(self):
        """Initialize the probabilities for replicon selection, given the abundances, etc.
        Normalization is always applied. Afterwards, produce a CDF which will be used for
        random sampling.
        """
        self.indexToName = {}
        prob = numpy.zeros(len(self.repliconRegistry))
        i = 0
        for repA in self.repliconRegistry.values():
            self.indexToName[i] = repA.name
            prob[i] = repA.parentCell.abundance / self.totalRawAbundance * repA.length()
            i += 1
        tp = sum(prob)
        self.pdf = numpy.divide(prob,tp)
        'Initialize the cumulative distribution function for the community replicons'
        self.cdf = numpy.hstack((0, numpy.cumsum(self.pdf)))
    
    def selectReplicon(self,x):
        'Return the index of a replicon by sampling CDF at given value x'
        idx = numpy.searchsorted(self.cdf, x) - 1
        if (idx < 0): #edge case when sample value is exactly zero.
            return 0
        return idx
    
    def pickReplicon(self,skipIndex=None):
        """Random selection of replicon from community, return the index"""
        if skipIndex == None:
            return self.selectReplicon(numpy.random.uniform())
        else:
            ri = skipIndex
            while (ri == skipIndex):
                ri = self.selectReplicon(numpy.random.uniform())
                print "picked and found",ri,"trying to avoid",skipIndex
            return ri
    
    def isSameReplicon(self):
        """Choose if the mate is intra or inter replicon associated. This is a simple
        binary paritioning with a chosen threshold frequency.
        """
        return numpy.random.uniform() < self.interRepliconProbability
    
    def isForwardStrand(self):
        'Coin flip for forward or reverse strand'
        return numpy.random.uniform() < 0.5
    
    def getRepliconByIndex(self, index):
        return self.repliconRegistry.get(self.indexToName[index])
    
    def getRepliconByName(self, name):
        return self.repliconRegistry.get(name)
    
    def unconstrainedReadLocation(self,repliconIndex):
        """Return a location (position and strand) on a replicon where the position is
        uniformly sampled across the replicons sequence.
        
        Returns tuple (pos=int, strand=bool)
        """
        replicon = self.getRepliconByIndex(repliconIndex)
        return (int(numpy.random.uniform() * replicon.length()), self.isForwardStrand())

    def constrainedReadLocation(self, repliconIndex, firstLocation, forward):
        """Return a location (position and strand) on a replicon where the position is
        constrained to follow a prescribed distribution relative to the position of the
        first location.
        
        Strand is always same as first.
        
        return location (pos=int, strand=bool)
        """
        replicon = self.getRepliconByIndex(repliconIndex)
        delta = drawDelta(500,replicon.length()-500)
        if forward: 
            # TODO The edge cases might have off-by-one errors, does it matter?'
            loc = firstLocation + delta
            if loc > replicon.length()-1:
                loc -= replicon.length()
        else:
            loc = firstLocation - delta
            if loc < 0:
                loc = replicon.length() - loc
        return (loc, forward)
    
    def getPDF(self):
        return self.pdf
    
    def getCDF(self):
        return self.cdf

def makeUnconstrainedPartA(otherRepIndex=None):
    repIdx = comm.pickReplicon(otherRepIndex)
    rep = comm.getRepliconByIndex(repIdx)
    pos1 = rep.randomCutSite('6cut')
    fwd = comm.isForwardStrand()
    pos2 = rep.nearestCutSiteAbove('4cut',pos1,fwd)
    seq = rep.subSeq(pos1,pos2,fwd)
    print rep,pos1,pos2,fwd
    return Part(seq, pos1, pos2, fwd, repIdx)

def makeConstrainedPartB(repIndex, firstLoc, forward):
    (loc, fwd) = comm.constrainedReadLocation(repIndex, firstLoc, forward)
    rep = comm.getRepliconByIndex(repIndex)
    pos1 = rep.nearestCutSiteByDistance('6cut',loc)
    pos2 = rep.nearestCutSiteAbove('4cut',pos1,fwd)
    seq = rep.subSeq(pos1,pos2,fwd)
    print rep,pos1,pos2,fwd
    return Part(seq,pos1,pos2,fwd,repIndex)

#
# Main 
#
if len(sys.argv) != 7:
    print('Usage: [n-frags] [read-len] [intra-prob] [table] [fasta] [output]')
    sys.exit(1)

# total number of reads to create
maxFragments = int(sys.argv[1])
# length of reads
readLength = int(sys.argv[2])

# Initialise community object
print "Initializing community"
comm = Community(float(sys.argv[3]),sys.argv[4],sys.argv[5])


hOutput = open(sys.argv[6],'wb')

fragCount = 0
while (fragCount < maxFragments):
    
    # Fragment creation
    
    # Create PartA
    # Steps
    # 1) pick a replicon at random
    # 2) pick a 6cutter site at random on replicon
    # 3) flip a coin for strand
    print "PARTA"
    partA = makeUnconstrainedPartA()
    
    # Create PartB
    # Steps
    # 1) choose if intra or inter replicon
    # 2) if intER create partB as above
    # 3) if intRA select from geometric
    
    print "PARTB"
    partB = None
    if comm.isSameReplicon():
        # intRA replicon
        print '\tINTRA todo'
        partB = makeConstrainedPartB(partA.repIdx, partA.pos1, partA.fwd)
        
    else:
        # intER replicon
        print '\tINTER'
        partB = makeUnconstrainedPartA(partA.repIdx)
        
    # Join parts
    # PartA + PartB
    fragment = partA.seq + partB.seq
    
    read1 = makeRead(fragment, True, readLength)
    read1.id = "frg" + str(fragCount) + "fwd"
    read1.description = partA.seq.id + " " + partA.seq.description

    read2 = makeRead(fragment, False, readLength)
    read2.id = "frg" + str(fragCount) + "rev"
    read2.description = partB.seq.id + " " + partB.seq.description
    
    writeReads(hOutput, [read1,read2], "fasta")
    
    fragCount += 1

hOutput.close()
