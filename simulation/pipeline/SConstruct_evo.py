from nestly import Nest
from nestly.scons import SConsWrap
import os
import os.path
import appconfig
import numpy as np

config = appconfig.read('config.yaml')

nest = Nest()
wrap = SConsWrap(nest, config['community']['folder'])
env = Environment(ENV=os.environ)

# Constants
wrap.add('seed', [config['seed']], create_dir=False)
wrap.add('genomes', [config['community']['seq']], create_dir=False)
wrap.add('basis_seq', [config['reference']['raw_seq']], create_dir=False)
wrap.add('seq_len', [config['reference']['seq_len']], create_dir=False)

# Variation
treeFolder = os.path.join(config['reference']['folder'], config['reference']['tree_folder'])
treePaths = appconfig.get_files(treeFolder, 'nwk')
wrap.add('tree', treePaths, label_func=os.path.basename)
wrap.add('branch_length', ['{:.4e}'.format(n) for n in np.logspace(-3, -6, num=10, endpoint=True).tolist()])

print treePaths

@wrap.add_target('generate_set')
def generate_set(outdir, c):
    tree = '{0[tree]}'.format(c)
    seq = '{0[basis_seq]}'.format(c)
    sources = [tree, seq]
    target = '{od}/{0[genomes]}'.format(c, od=outdir)
    action = 'bin/pbsrun_SGEVOLVER.sh ' \
             '{0[seed]} {0[branch_length]} {0[seq_len]} $SOURCES.abspath $TARGET.abspath'.format(c)
    return 'hr', env.Command(target, sources, action)

wrap.add_controls(Environment())
