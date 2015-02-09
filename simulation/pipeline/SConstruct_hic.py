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
commPaths = appconfig.get_communities(config)
wrap.add('community', commPaths, label_func=os.path.basename)


index_suffixes = ['.amb', '.ann', '.bwt', '.pac', '.sa']
@wrap.add_target('index_ref')
def index_ref(outdir, c):
    source = '{0[community]}/{0[refseq]}'.format(c)
    target = [source + suf for suf in index_suffixes]
    action = 'bin/pbsrun_INDEX.sh $SOURCE.abspath'.format(c)
    return env.Command(target, source, action)

wrap.add('comm_table', [config['community']['table']], label_func=stripext)
wrap.add('hic_n_frag', config['hic_n_frag'])


@wrap.add_target('generate_hic')
@name_targets
def generate_hic(outdir, c):
    source = '{0[community]}/{0[refseq]}'.format(c)
    target = '{od}/{0[hic_base]}.fasta'.format(c, od=outdir)
    action = 'bin/pbsrun_HiC.sh ' \
             '{0[seed]} {0[hic_n_frag]} {0[hic_read_length]} {0[hic_inter_prob]} ' \
             '{0[community]}/{0[comm_table]} $SOURCE.abspath $TARGET.abspath'.format(c)
    return 'hr', env.Command(target, source, action)

wrap.add_controls(Environment())

