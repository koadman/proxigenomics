from nestly import Nest, stripext
from nestly.scons import SConsWrap, name_targets
import os
import os.path
import appconfig

config = appconfig.read('config.yaml')

nest = Nest()
wrap = SConsWrap(nest, config['hic_folder'])
env = Environment(ENV=os.environ)

# Constants
wrap.add('seed', [config['seed']], create_dir=False)
wrap.add('refseq', [config['community']['seq']], create_dir=False)
wrap.add('hic_base', [config['hic_base']], create_dir=False)
wrap.add('hic_inter_prob', [config['hic_inter_prob']], create_dir=False)
wrap.add('hic_read_length', [config['hic_read_length']], create_dir=False)

# Variation
genomes = appconfig.find_files(config['community']['folder'], config['community']['seq'])
commPaths = [os.path.dirname(pn) for pn in genomes]
wrap.add('community', commPaths)

tableFolder = os.path.join(config['reference']['folder'], config['reference']['table_folder'])
tablePaths = appconfig.get_files(tableFolder, 'table')
wrap.add('comm_table', tablePaths, label_func=os.path.basename)

wrap.add('hic_n_frag', config['hic_n_frag'])


@wrap.add_target('generate_hic')
def generate_hic(outdir, c):
    source = '{1[community][folder]}/{0[community]}/{0[refseq]}'.format(c, config)
    target = '{od}/{0[hic_base]}.fasta'.format(c, od=outdir)
    action = 'bin/pbsrun_HiC.sh ' \
             '{0[seed]} {0[hic_n_frag]} {0[hic_read_length]} {0[hic_inter_prob]} ' \
             '{0[comm_table]} $SOURCE.abspath $TARGET.abspath'.format(c)
    return env.Command(target, source, action)

wrap.add_controls(Environment())

