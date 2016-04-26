from nestly import Nest, stripext
from nestly.scons import SConsWrap, name_targets
import os
import os.path
import appconfig

config = appconfig.read('config.yaml')

nest = Nest()
wrap = SConsWrap(nest, config['map_folder'])
env = Environment(ENV=os.environ)

# Used for resolving what type of execution environment will be used.
exec_env = appconfig.ExecutionEnvironment(ARGUMENTS, supported_env=['pbs', 'sge'])

# Constant
wrap.add('refseq', [config['community']['seq']], create_dir=False)

# Variation
genomes = appconfig.find_files(config['community']['folder'], config['community']['seq'])
# For testing - restrict to only star topology
#commPaths = [os.path.dirname(pn) for pn in genomes if 'ladder' not in pn]
commPaths = [os.path.dirname(pn) for pn in genomes]
wrap.add('community', commPaths)

tableFolder = os.path.join(config['reference']['folder'], config['reference']['table_folder'])
tablePaths = appconfig.get_files(tableFolder, 'table')
wrap.add('comm_table', tablePaths, label_func=os.path.basename)

wrap.add('wgs_xfold', config['wgs_xfold'])


@wrap.add_target('make_ctg2ref')
@name_targets
def make_ctg2ref(outdir, c):
    com = c['community']
    tbl = os.path.basename(c['comm_table'])

    query = os.path.join(
                os.path.abspath(config['wgs_folder']), com, tbl,
                str(c['wgs_xfold']), config['wgs_asmdir'],
                '{0[wgs_base]}.contigs.fasta'.format(config))

    subject = os.path.join(config['community']['folder'], com, config['community']['seq'])

    source = [subject, query]
    target = os.path.join(outdir, config['ctg2ref'])

    action = exec_env.resolve_action({
        'pbs': 'bin/pbsrun_LAST.sh $SOURCES.abspath $TARGET.abspath',
        'sge': 'bin/sgerun_LAST.sh $SOURCES.abspath $TARGET.abspath'
    })

    return 'output', env.Command(target, source, action)

@wrap.add_target('make_truth')
@name_targets
def make_truth(outdir, c):
    source = str(c['make_ctg2ref']['output'])
    target = os.path.join(outdir, config['truth_table'])

    action = exec_env.resolve_action({

        'pbs': 'bin/pbsrun_MKTRUTH.sh '
               '{1[ctg_afmt]} {1[ctg_ofmt]} {1[ctg_minlen]} {1[ctg_mincov]} {1[ctg_minid]} '
               '$SOURCES.abspath $TARGET.abspath'.format(c, config),

        'sge': 'bin/sgerun_MKTRUTH.sh '
               '{1[ctg_afmt]} {1[ctg_ofmt]} {1[ctg_minlen]} {1[ctg_mincov]} {1[ctg_minid]} '
               '$SOURCES.abspath $TARGET.abspath'.format(c, config)
    })

    return 'output', env.Command(target, source, action)


@wrap.add_target('make_wgs2ctg')
@name_targets
def make_wgs2ctg(outdir, c):
    com = c['community']
    tbl = os.path.basename(c['comm_table'])

    # TODO find a better way to obtain the path to WGS reads
    query = appconfig.get_wgs_reads(
                os.path.join(os.path.abspath(config['wgs_folder']),
                com, tbl, str(c['wgs_xfold'])),
                config)

    # TODO likewise
    subject = os.path.join(
                os.path.abspath(config['wgs_folder']),
                com, tbl, str(c['wgs_xfold']), config['wgs_asmdir'],
                '{0[wgs_base]}.contigs.fasta'.format(config))

    target = os.path.join(outdir, config['wgs2ctg'])
    source = [subject] + query

    action = exec_env.resolve_action({
        'pbs': 'bin/pbsrun_MAP.sh $SOURCES.abspath $TARGET.abspath',
        'sge': 'bin/sgerun_MAP.sh $SOURCES.abspath $TARGET.abspath'
    })

    return 'output', env.Command(target, source, action)


wrap.add('hic_n_frag', config['hic_n_frag'])

@wrap.add_target('make_hic2ctg')
@name_targets
def make_hic2ctg(outdir, c):
    com = c['community']
    tbl = os.path.basename(c['comm_table'])

    query = os.path.join(os.path.abspath(config['hic_folder']), com, tbl,
                         str(c['hic_n_frag']), '{0[hic_base]}.fasta'.format(config))
    subject = os.path.join(
                os.path.abspath(config['wgs_folder']), com, tbl,
                str(c['wgs_xfold']), config['wgs_asmdir'],
                '{0[wgs_base]}.contigs.fasta'.format(config))
    source = [subject, query]
    target = os.path.join(outdir, config['hic2ctg'])

    action = exec_env.resolve_action({
        'pbs': 'bin/pbsrun_MAP.sh $SOURCES.abspath $TARGET.abspath',
        'sge': 'bin/sgerun_MAP.sh $SOURCES.abspath $TARGET.abspath'
    })

    return 'output', env.Command(target, source, action)


wrap.add('hic_min_cov', [0, 0.5, 0.95])
wrap.add('hic_min_qual', [0, 30, 60])


@wrap.add_target('filter_hic2ctg')
@name_targets
def filter_hic(outdir, c):
    source = str(c['make_hic2ctg']['output'])
    target = os.path.join(outdir, config['hic2ctg'])
    action = exec_env.resolve_action({
        'pbs': 'bin/pbsrun_BAMFILTER.sh {0[hic_min_cov]} {0[hic_min_qual]} $SOURCE $TARGET'.format(c),
        'sge': 'bin/sgerun_BAMFILTER.sh {0[hic_min_cov]} {0[hic_min_qual]} $SOURCE $TARGET'.format(c)
    })
    return 'output', env.Command(target, source, action)


wrap.add_controls(Environment())
