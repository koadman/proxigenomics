from nestly import Nest, stripext
from nestly.scons import SConsWrap, name_targets
import os
import os.path
import appconfig
import numpy as np

config = appconfig.read('config.yaml')

nest = Nest()
wrap = SConsWrap(nest, config['evo_folder'])
env = Environment(ENV=os.environ)

# Constants
wrap.add('seed', [config['seed']], create_dir=False)
wrap.add('genomes', [config['community']['seq']], create_dir=False)
wrap.add('basis_seq', [config['community']['basis']], create_dir=False)
wrap.add('seq_len', [config['community']['seq_len']], create_dir=False)
# Variation
wrap.add('tree', config['community']['tree'], label_func=stripext)
wrap.add('comm_table', [config['community']['table']], label_func=stripext)
wrap.add('branch_length', np.logspace(-3, -6, num=10, endpoint=True).tolist())


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

