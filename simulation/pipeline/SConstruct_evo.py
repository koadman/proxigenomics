from nestly import Nest
from nestly.scons import SConsWrap
import os
import os.path
import appconfig
import numpy as np
import math

config = appconfig.read('config.yaml')

nest = Nest()
wrap = SConsWrap(nest, config['community']['folder'])
env = Environment(ENV=os.environ)

# Constants
wrap.add('seed', [config['seed']], create_dir=False)
wrap.add('genomes', [config['community']['seq']], create_dir=False)
wrap.add('basis_seq', [config['reference']['raw_seq']], create_dir=False)
wrap.add('seq_len', [config['reference']['seq_len']], create_dir=False)
wrap.add('sg_scale', [config['reference']['sg_scale']], create_dir=False)

# Variation

# trees are sourced from folder
seqFolder = os.path.join(config['reference']['folder'])
treeFolder = os.path.join(config['reference']['folder'], config['reference']['tree_folder'])
treePaths = appconfig.get_files(treeFolder, config['reference']['tree_suffix'])
wrap.add('tree', treePaths, label_func=os.path.basename)

# scale divergence of tree evenly across log space
wrap.add('branch_length', ['{:.4e}'.format(n) for n in np.logspace(
    start=math.log(config['reference']['tree_scale']['max']),
    stop=math.log(config['reference']['tree_scale']['min']),
    num=config['reference']['tree_scale']['steps'],
    base=math.e,
    endpoint=True).tolist()])


@wrap.add_target('generate_set')
def generate_set(outdir, c):
    tree = os.path.join(treeFolder, '{0[tree]}'.format(c))
    seq = os.path.join(seqFolder, '{0[basis_seq]}'.format(c))
    sources = [tree, seq]
    target = '{od}/{0[genomes]}'.format(c, od=outdir)
    action = 'bin/pbsrun_SGEVOLVER.sh ' \
             '{0[seed]} {0[branch_length]} {0[sg_scale]} {0[seq_len]} $SOURCES.abspath $TARGET.abspath'.format(c)
    return 'hr', env.Command(target, sources, action)


@wrap.add_target('reconstruct')
def reconstruct(outdir, c):
    base = 'reconstruct'
    source = '{od}/{0[genomes]}'.format(c, od=outdir)
    target = '{od}/{base}/ani/perc_ids.tab'.format(od=outdir,base=base)
    action = 'bin/pbsrun_MKTREE.sh $SOURCE.abspath {od}/{base}'.format(c, od=outdir, base=base)
    return 'recon', env.Command(target, source, action)

wrap.add_controls(Environment())
