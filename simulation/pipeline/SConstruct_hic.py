from nestly import Nest, stripext
from nestly.scons import SConsWrap, name_targets
import glob
import os
import os.path

nest = Nest()
wrap = SConsWrap(nest, 'hic_data')
env = Environment(ENV = os.environ)

commPaths = [os.path.abspath(dn) for dn in glob.glob('references/*')]

# Constants
wrap.add('seed', [2136841], create_dir=False)
wrap.add('refseq', ['genomes.fasta'], create_dir=False)
wrap.add('hic_base', ['hic'], create_dir=False)
wrap.add('hic_inter_prob', [0.9], create_dir=False)
wrap.add('hic_read_length', [150], create_dir=False)

# Variation
wrap.add('community', commPaths, label_func=os.path.basename)

index_suffixes = ['.amb', '.ann', '.bwt', '.pac', '.sa']
@wrap.add_target('index_ref')
def index_ref(outdir, c):
	source = '{0[community]}/{0[refseq]}'.format(c)
	target = [source + suf for suf in index_suffixes]
	action = 'cd {od} && /panfs/panspermia/120274/work/hi-c/sim/bin/pbsrun_INDEX.sh $SOURCE.abspath'.format(c,od=outdir)
	return env.Command(target,source,action)

wrap.add('comm_table', ['uniform.table'], label_func=stripext)
wrap.add('hic_n_frag', [100000,250000,500000])

@wrap.add_target('generate_hic')
@name_targets
def generate_hic(outdir, c):
	target = '{od}/{0[hic_base]}.fasta'.format(c,od=outdir)
	source = '{0[community]}/{0[refseq]}'.format(c)
	action = 'cd {od} && /panfs/panspermia/120274/work/hi-c/sim/bin/pbsrun_HiC.sh {0[seed]} {0[hic_n_frag]} {0[hic_read_length]} {0[hic_inter_prob]} {0[community]}/{0[comm_table]} $SOURCE.abspath {0[hic_base]}'.format(c,od=outdir)
	return 'hr',env.Command(target,source,action)

wrap.add_controls(Environment())

