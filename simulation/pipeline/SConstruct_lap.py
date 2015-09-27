from nestly import Nest, stripext
from nestly.scons import SConsWrap, name_targets
import os
import os.path
import appconfig

config = appconfig.read('config.yaml')

nest = Nest()
wrap = SConsWrap(nest, config['lap_folder'])
env = Environment(ENV=os.environ)

# Used for resolving what type of execution environment will be used.
exec_env = appconfig.ExecutionEnvironment(ARGUMENTS, supported_env=['pbs', 'sge'])

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

    src_dir = os.path.join(config['wgs_folder'], '/'.join(outdir.split('/')[1:]))
    sources = ['{1[community][folder]}/{0[community]}/{0[refseq]}'.format(c, config)] + \
              appconfig.get_wgs_reads(src_dir, config)

    target = os.path.join(outdir, 'lap_ref.prob')

    action = exec_env.resolve_action({
        'pbs': 'bin/pbsrun_LAP.sh $SOURCES.abspath $TARGET.abspath'.format(c, od=outdir),
        'sge': 'bin/sgerun_LAP.sh $SOURCES.abspath $TARGET.abspath'.format(c, od=outdir)
    })

    return 'lap_ref', env.Command(target, sources, action)

@wrap.add_target('reads2ctg')
@name_targets
def generate_wgs(outdir, c):

    src_dir = os.path.join(config['wgs_folder'], '/'.join(outdir.split('/')[1:]))
    sources = ['{1}/{2}/{0[wgs_base]}.contigs.fasta'.format(c, src_dir, config['wgs_asmdir'])] + \
              appconfig.get_wgs_reads(src_dir, config)

    target = os.path.join(outdir, 'lap_ctg.prob')

    action = exec_env.resolve_action({
        'pbs': 'bin/pbsrun_LAP.sh $SOURCES.abspath $TARGET.abspath'.format(c, od=outdir),
        'sge': 'bin/sgerun_LAP.sh $SOURCES.abspath $TARGET.abspath'.format(c, od=outdir)
    })

    return 'lap_ctg', env.Command(target, sources, action)

@wrap.add_target('sumprob')
@name_targets
def sumprob(outdir, c):

    sources = [str(c['reads2ref']['lap_ref']), str(c['reads2ctg']['lap_ctg'])]
    targets = [os.path.join(outdir, 'lap_ref.sum'), os.path.join(outdir, 'lap_ctg.sum')]

    action = exec_env.resolve_action({
        'pbs': 'bin/pbsrun_LAPSUM.sh $SOURCES.abspath $TARGETS.abspath'.format(c, od=outdir),
        'sge': 'bin/sgerun_LAPSUM.sh $SOURCES.abspath $TARGETS.abspath'.format(c, od=outdir)
    })

    return 'sum_ref', 'sub_ctg', env.Command(targets, sources, action)

wrap.add_controls(Environment())
