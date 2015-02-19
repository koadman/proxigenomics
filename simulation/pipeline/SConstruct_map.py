from nestly import Nest, stripext
from nestly.scons import SConsWrap, name_targets
import os
import os.path
import numpy
import appconfig

#
# Helper functions
#


def pick_alignment(fname, c):
    return [fn for fn in c['align_files'] if fn.endswith(fname)]


def get_wgs_fasta(c):
    com = os.path.basename(c['community'])
    tbl = stripext(c['hic_table'])
    return os.path.join(
            os.path.abspath(config['wgs_folder']), com, tbl,
            str(c['wgs_xfold']), config['wgs_asmdir'],
            '{0}.contigs.fasta'.format(config['wgs_base']))

import types
def prepend_paths(path, fnames):
    if isinstance(fnames, types.StringTypes):
        fnames = [fnames]
    return [os.path.join(path, fn) for fn in fnames]


config = appconfig.read('config.yaml')


nest = Nest()
wrap = SConsWrap(nest, config['map_folder'])
env = Environment(ENV=os.environ)

# Constant
wrap.add('refseq', [config['community']['seq']], create_dir=False)

commPaths = appconfig.get_communities(config)
wrap.add('community', commPaths, label_func=os.path.basename)

wrap.add('com_table', [config['community']['table']], label_func=stripext)
wrap.add('hic_n_frag', config['hic_n_frag'])
wrap.add('wgs_xfold', config['wgs_xfold'])

wrap.add_aggregate('align_files', list)
wrap.add_aggregate('graph_files', list)
wrap.add_aggregate('cl_input', list)


@wrap.add_target('make_hic2ctg')
def map_hic2ctg(outdir, c):
    com = os.path.basename(c['community'])
    tbl = stripext(c['com_table'])

    query = os.path.join(os.path.abspath(config['hic_folder']), com, tbl,
                         str(c['hic_n_frag']), '{0[hic_base]}.fasta'.format(config))

    subject = os.path.join(
                os.path.abspath(config['wgs_folder']), com, tbl,
                str(c['wgs_xfold']), config['wgs_asmdir'],
                '{0[wgs_base]}.contigs.fasta'.format(config))

    target = os.path.join(outdir, config['hic2ctg'])

    source = [subject, query]

    action = 'bin/pbsrun_MAP.sh $SOURCES.abspath $TARGET.abspath'

    c['align_files'].append(target)

    return env.Command(target, source, action)


@wrap.add_target('make_ctg2ref')
@name_targets
def map_ctg2ref(outdir, c):
    com = os.path.basename(c['community'])
    tbl = stripext(c['com_table'])

    query = os.path.join(
                os.path.abspath(config['wgs_folder']), com, tbl,
                str(c['wgs_xfold']), config['wgs_asmdir'],
                '{0[wgs_base]}.contigs.fasta'.format(config))

    subject = os.path.join(c['community'], config['community']['seq'])

    target = os.path.join(outdir, config['ctg2ref'])

    source = [subject, query]

    action = 'bin/pbsrun_LAST.sh $SOURCES.abspath $TARGET.abspath'

    return 'output', env.Command(target, source, action)


@wrap.add_target('make_wgs2ctg')
def map_wgs2ctg(outdir, c):
    com = os.path.basename(c['community'])
    tbl = stripext(c['com_table'])

    # TODO find a better way to obtain the path to WGS reads
    query = appconfig.get_wgs_reads(
                os.path.join(os.path.abspath(config['wgs_folder']), com, tbl, str(c['wgs_xfold'])),
                config)

    # TODO likewise
    subject = os.path.join(
                os.path.abspath(config['wgs_folder']),
                com, tbl, str(c['wgs_xfold']), config['wgs_asmdir'],
                '{0[wgs_base]}.contigs.fasta'.format(config))

    target = os.path.join(outdir, config['wgs2ctg'])
    source = [subject] + query
    action = 'bin/pbsrun_MAP.sh $SOURCES.abspath $TARGET.abspath'

    c['align_files'].append(target)

    return env.Command(target, source, action)


@wrap.add_target('make_truth')
@name_targets
def make_truth(outdir, c):
    source = str(c['make_ctg2ref']['output'])
    target = os.path.join(outdir, config['truth_table'])
    action = 'bin/alignmentToTruth.py ' \
             '--afmt {1[ctg_afmt]} ' \
             '--ofmt {1[ctg_ofmt]} ' \
             '--minlen {1[ctg_minlen]} ' \
             '--mincov {1[ctg_mincov]} ' \
             '--minid {1[ctg_minid]} ' \
             '$SOURCES.abspath $TARGET.abspath'.format(c, config)
    return 'truth', env.Command(target, source, action)


@wrap.add_target('make_graph')
def make_graph(outdir, c):
    hic_sam = pick_alignment(config['hic2ctg'], c)
    wgs_bam = [os.path.splitext(fn)[0] + ".bam" for fn in c['align_files'] if fn.endswith(config['wgs2ctg'])]
    source = hic_sam + wgs_bam
    target = prepend_paths(outdir, ['edges.csv', 'nodes.csv'])
    c['graph_files'].append(target)
    action = 'bin/pbsrun_GRAPH.sh $SOURCES.abspath $TARGETS.abspath'
    return env.Command(target, source, action)

#
#  Everything below here is MCL specific but should be made agnostic of algorithm
#  or deal with multiple algorithms in "cluster_method".
#
wrap.add('cluster_method', ['mcl'])

@wrap.add_target('make_cluster_input')
def make_cluster_input(outdir, c):
    source = c['graph_files']
    target = prepend_paths(outdir, config['cluster']['input'])
    action = 'bin/pbsrun_MKMCL.sh {1[ctg_minlen]} $SOURCES.abspath $TARGET.abspath'.format(c, config)
    c['cl_input'].append(target)
    return env.Command(target, source, action)


wrap.add('mcl_inflation', numpy.linspace(1.1, 2.0, 9))

@wrap.add_target('do_mcl')
def do_mcl(outdir, c):
    # TODO run over both weighted/unweighted?
    source = c['cl_input'][0][0]
    target = prepend_paths(outdir, config['cluster']['output'])
    action = 'bin/pbsrun_MCL.sh {0[mcl_inflation]} $SOURCE.abspath $TARGET.abspath'.format(c)
    return env.Command(target, source, action)


@wrap.add_target('do_score')
def do_score(outdir, c):
    cl_out = os.path.join(outdir, config['cluster']['output'])
    ttable = str(c['make_truth']['truth'])
    # this target consumes truth table and clustering output
    source = [ttable, cl_out]
    # this target creates 3 output files
    target = ['{0}.{1}'.format(cl_out, suffix) for suffix in ['joined', 'f1', 'vm']]
    action = 'bin/pbsrun_SCORE.sh $SOURCES.abspath'
    return env.Command(target, source, action)


wrap.add_controls(Environment())
