from nestly import Nest, stripext
from nestly.scons import SConsWrap, name_targets
import os
import os.path
import appconfig

config = appconfig.read('config.yaml')

nest = Nest()
wrap = SConsWrap(nest, config['quast_folder'])
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
wrap.add('community', commPaths)

tableFolder = os.path.join(config['reference']['folder'], config['reference']['table_folder'])
tablePaths = appconfig.get_files(tableFolder, 'table')
wrap.add('comm_table', tablePaths, label_func=os.path.basename)

wrap.add('wgs_xfold', config['wgs_xfold'])

@wrap.add_target('run_quast')
@name_targets
def generate_wgs(outdir, c):

    asm_path = os.path.join(config['wgs_folder'], '/'.join(outdir.split('/')[1:]), config['wgs_asmdir'])
    ctg_file = '{0}/{1}.contigs.fasta'.format(asm_path, config['wgs_base'])
    sources = ['{1[community][folder]}/{0[community]}/{0[refseq]}'.format(c, config), ctg_file]
    target = os.path.join(outdir, 'quast/report.tsv')

    action = exec_env.resolve_action({
        'sge': 'bin/sgerun_QUAST.sh $SOURCES.abspath $TARGET.abspath'.format(c, od=outdir)
    })

    return 'report', env.Command(target, sources, action)

wrap.add_controls(Environment())
