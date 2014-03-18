from Bio import SeqIO
import re
import sys

matcher = re.compile(r'([0-9]+)M')
def countAligned(cigar):
    return sum([int(s) for s in matcher.findall(fields[5])])

class Alignment:
    def __init__(self,bases=0):
        self.bases =  bases
    def __repr__(self):
        return repr((self.bases))
    def __str__(self):
        return str(self.bases)
    def addbases(self,bases):
        self.bases += bases

if len(sys.argv) != 4:
	print 'Usage: [min length] [query fasta] [sam file]'
	sys.exit(1)

minLength = int(sys.argv[1])
seqLength = {rec.id: len(rec) for rec in SeqIO.parse(sys.argv[2],'fasta')}

taxons = {}
i=0
hin = open(sys.argv[3],'r')
for l in hin:
    if l.startswith('@'): continue
    fields = l.rsplit()
    tx = taxons.get(fields[2])
    if tx is None:
        tx = {}
        taxons[fields[2]] = tx
    txal = tx.get(fields[0])
    if txal is None:
        txal = Alignment()
        tx[fields[0]] = txal
    txal.addbases(countAligned(fields[5]))
hin.close()

def tostring(d):
    s = ''
    for k,v in d.iteritems():
        s += '{k} {v} '.format(k=str(k), v=str(v))
    return s

#
# We need to decide on assignment.
#
# There are contigs which are assigned to more than one taxon.
# What shall be done with these? Simply taking the largest alignmetn
# as the winner might be acceptable for large differences, but I've
# seen examples where both are significant alignments and picking one
# is quite arbitrary.
#

# PICK THE WINNER
# For each scaffold, determine the best alignment by length. The alignment subject
# will then be considered the true source.
best = {}
for t,s in taxons.iteritems():
    for sn,aln in s.iteritems():
		bi = best.get(sn)
		if bi is None or aln.bases > bi['aln'].bases:
			best[sn] = {'tx': t, 'aln': aln, 'slen': seqLength[sn]}

# Write out the list of winners, with additional columns for validation
for k,v in best.iteritems():
	if v['aln'].bases > minLength:
		print '{seq_name} {tax_name} {cov:.4}'.format(
			seq_name=k, tax_name=v['tx'], cov=float(v['aln'].bases)/float(v['slen']))
	