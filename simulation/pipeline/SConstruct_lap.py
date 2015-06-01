from nestly import Nest, stripext
from nestly.scons import SConsWrap, name_targets
import os
import os.path
import appconfig

config = appconfig.read('config.yaml')

nest = Nest()
wrap = SConsWrap(nest, config['lap_folder'])
env = Environment(ENV=os.environ)

# Constants
wrap.add('seed', [config['seed']], create_dir=False)
wrap.add('refseq', [config['community']['seq']], create_dir=False)
wrap.add('wgs_read_length', [config['wgs_read_length']], create_dir=False)
wrap.add('wgs_insert_length', [config['wgs_insert_length']], create_dir=False)
wrap.add('wgs_insert_sd', [config['wgs_insert_sd']], create_dir=False)
wrap.add('wgs_base', [config['wgs_base']],  create_dir=False)

# Variation
genomes = appconfig.find_files(config['community']['folder'], config['community']['seq'])
commPaths = [os.path.dirname(pn) for pn in genomes]
# For testing - constraint on branch length
#commPaths = [pn for pn in commPaths if float(os.path.basename(pn)) > 0.4 and float(os.path.basename(pn)) < 1]
wrap.add('community', commPaths)

tableFolder = os.path.join(config['reference']['folder'], config['reference']['table_folder'])
tablePaths = appconfig.get_files(tableFolder, 'table')
wrap.add('comm_table', tablePaths, label_func=os.path.basename)

wrap.add('wgs_xfold', config['wgs_xfold'])

@wrap.add_target('reads2ref')
@name_targets
def generate_wgs(outdir, c):

    sources = ['{1[community][folder]}/{0[community]}/{0[refseq]}'.format(c, config)] + \
              appconfig.get_wgs_reads(config['wgs_folder'], config)

    target = os.path.join(outdir, 'lap_ref.prob')

    action = 'bin/pbsrun_LAP.sh $SOURCES.abspath $TARGET.abspath'.format(c, od=outdir)
    print sources, target
    return 'lap_ref', env.Command(target, sources, action)

@wrap.add_target('reads2ctg')
@name_targets
def generate_wgs(outdir, c):

    sources = ['{1[wgs_folder]}/{0[wgs_base]}.contigs.fasta'.format(c, config)] + \
              appconfig.get_wgs_reads(config['wgs_folder'], config)

    target = os.path.join(outdir, 'lap_ctg.prob')

    action = 'bin/pbsrun_LAP.sh $SOURCES.abspath $TARGET.abspath'.format(c, od=outdir)
    print sources, target
    return 'lap_ctg', env.Command(target, sources, action)

import sys
sys.exit(1)

wrap.add_controls(Environment())
