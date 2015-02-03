from nestly import Nest, stripext
from nestly.scons import SConsWrap, name_targets
import glob
import os
import os.path

nest = Nest()
wrap = SConsWrap(nest, 'wgs_data')
env = Environment(ENV=os.environ)

commPaths = [os.path.abspath(dn) for dn in glob.glob('references/*')]

# Variation
wrap.add('community', commPaths, label_func=os.path.basename)
wrap.add('comm_table', ['uniform.table'], label_func=stripext)
wrap.add('wgs_xfold', [1, 2, 5])

# Constants
wrap.add('seed', [2136841], create_dir=False)
wrap.add('refseq', ['genomes.fasta'], create_dir=False)
wrap.add('wgs_read_length', [150], create_dir=False)
wrap.add('wgs_insert_length', [450], create_dir=False)
wrap.add('wgs_insert_sd', [100], create_dir=False)
wrap.add('wgs_base', ['wgs'],  create_dir=False)


@wrap.add_target('generate_wgs')
@name_targets
def generate_wgs(outdir, c):
    target = [os.path.join(outdir, of) for of in ['wgs1.fq', 'wgs2.fq']]
    source = '{0[community]}/{0[refseq]}'.format(c)
    action = 'bin/pbsrun_ART.sh {0[seed]} {0[wgs_insert_length]} {0[wgs_insert_sd]} ' \
             '{0[wgs_read_length]} {0[wgs_xfold]} $SOURCE.abspath {od}/{0[wgs_base]}'.format(c, od=outdir)
    return 'r1', 'r2', env.Command(target, source, action)


@wrap.add_target('assemble_wgs')
@name_targets
def assemble_wgs(work_dir, c):
    asm_dir = work_dir + '/asm'.format(c)
    target = '{od}/{0[wgs_base]}.contigs.fasta'.format(c, od=asm_dir)
    source = [str(c['generate_wgs']['r1']), str(c['generate_wgs']['r2'])]
    action = 'bin/a5submit.sh -mwf -t {0[wgs_base]} $SOURCES.abspath {od}'.format(c, od=asm_dir)
    return 'ctg', env.Command(target, source, action)


index_suffixes = ['.amb', '.ann', '.bwt', '.pac', '.sa']
@wrap.add_target('index_ctg')
def index_ctg(outdir, c):
    source = str(c['assemble_wgs']['ctg'])
    target = [source + suf for suf in index_suffixes]
    action = 'bin/pbsrun_INDEX.sh $SOURCE.abspath'.format(c)
    return env.Command(target, source, action)

wrap.add_controls(Environment())

