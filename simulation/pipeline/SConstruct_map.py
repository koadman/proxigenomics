from nestly import Nest, stripext
from nestly.scons import SConsWrap, name_targets
import fnmatch
import glob
import os
import os.path
import numpy
import appconfig

#
# Helper functions
#


def pick_sam(fname, c):
    return [fn for fn in c['sam_files'] if fn.endswith(fname)]


def get_wgs_fasta(c):
    com = os.path.basename(c['community'])
    tbl = stripext(c['hic_table'])
    return os.path.join(
            os.path.abspath(config['wgs_folder']), com, tbl,
            str(c['wgs_xfold']), config['wgs_asmdir'],
            '{0}.contigs.fasta'.format(config['wgs_base']))


def prepend_paths(path, fnames):
    return [os.path.join(path, fn) for fn in fnames]


config = appconfig.read('config.yaml')


nest = Nest()
wrap = SConsWrap(nest, config['map_folder'])
env = Environment(ENV=os.environ)

commPaths = appconfig.get_communities(config)
wrap.add('community', commPaths, label_func=os.path.basename)

wrap.add('hic_table', [config['community']['table']], label_func=stripext)
wrap.add('hic_n_frag', config['hic_n_frag'])
wrap.add('wgs_xfold', config['wgs_xfold'])

wrap.add_aggregate('sam_files', list)
wrap.add_aggregate('graph_files', list)
wrap.add_aggregate('cl_input', list)
wrap.add_aggregate('truth', list)


@wrap.add_target('make_hic2ctg')
def map_hic2ctg(outdir, c):
    com = os.path.basename(c['community'])
    tbl = stripext(c['hic_table'])

    query = os.path.join(os.path.abspath(config['hic_folder']), com, tbl,
                         str(c['hic_n_frag']), '{0[hic_base]}.fasta'.format(config))

    subject = os.path.join(
                os.path.abspath(config['wgs_folder']), com, tbl,
                str(c['wgs_xfold']), config['wgs_asmdir'],
                '{0[wgs_base]}.contigs.fasta'.format(config))

    target = os.path.join(outdir, config['hic2ctg'])

    source = [subject, query]

    action = 'bin/pbsrun_MAP.sh $SOURCES.abspath $TARGET.abspath'

    c['sam_files'].append(target)

    return env.Command(target, source, action)


@wrap.add_target('make_ctg2ref')
def map_ctg2ref(outdir, c):
    com = os.path.basename(c['community'])
    tbl = stripext(c['hic_table'])

    query = os.path.join(
                os.path.abspath(config['wgs_folder']), com, tbl,
                str(c['wgs_xfold']), config['wgs_asmdir'],
                '{0[wgs_base]}.contigs.fasta'.format(config))

    subject = os.path.join(c['community'], config['community']['seq'])

    target = os.path.join(outdir, config['ctg2ref'])

    source = [subject, query]

    action = 'bin/pbsrun_MAP.sh $SOURCES.abspath $TARGET.abspath'

    c['sam_files'].append(target)

    return env.Command(target, source, action)


@wrap.add_target('make_wgs2ctg')
def map_wgs2ctg(outdir, c):
    com = os.path.basename(c['community'])
    tbl = stripext(c['hic_table'])

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
    action = 'bin/pbsrun_MAP2.sh $SOURCES.abspath $TARGET.abspath'

    c['sam_files'].append(target)

    return env.Command(target, source, action)


@wrap.add_target('make_sam2bam')
def make_bam(outdir, c):
    target = [os.path.splitext(dn)[0] + ".bam" for dn in c['sam_files']]
    source = c['sam_files']
    action = 'bin/pbsrun_SAMTOOLS.sh {0}'.format(outdir)
    return env.Command(target, source, action)


@wrap.add_target('make_graph')
def make_graph(outdir, c):
    hic_sam = pick_sam(config['hic2ctg'], c)
    wgs_bam = [os.path.splitext(fn)[0] + ".bam" for fn in c['sam_files'] if fn.endswith(config['wgs2ctg'])]
    source = hic_sam + wgs_bam
    target = prepend_paths(outdir, ['edges.csv', 'nodes.csv'])
    c['graph_files'].append(target)
    action = 'bin/pbsrun_GRAPH.sh $SOURCES.abspath $TARGETS.abspath'
    return env.Command(target, source, action)


wrap.add('cluster_method', ['mcl'])
wrap.add('score_min_length', [500, 1000])


@wrap.add_target('make_truth')
def make_truth(outdir, c):
    seq = get_wgs_fasta(c)
    source = [seq] + pick_sam(config['ctg2ref'], c)
    target = os.path.join(outdir, config['truth_table'])
    c['truth'].append(target)
    action = 'bin/parseSamCigar.py {0[score_min_length]} $SOURCES.abspath $TARGET.abspath'.format(c)
    return env.Command(target, source, action)


@wrap.add_target('make_cluster_input')
def make_cluster_input(outdir, c):
    source = c['graph_files']
    target = prepend_paths(outdir, ['mclIn.weighted'])
    action = 'bin/makeMCLinput.py {0[score_min_length]} $SOURCES.abspath $TARGET.abspath'.format(c)
    c['cl_input'].append(target)
    return env.Command(target, source, action)


wrap.add('mcl_inflation', numpy.linspace(1.1, 2.0, 9))


@wrap.add_target('do_mcl')
def do_mcl(outdir, c):
    # TODO run over both weighted/unweighted
    source = c['cl_input'][0][0]
    target = prepend_paths(outdir, ['cl_out'])
    action = 'bin/pbsrun_MCL.sh {0[mcl_inflation]} $SOURCE.abspath $TARGET.abspath'.format(c)
    return env.Command(target, source, action)


@wrap.add_target('do_score')
def do_score(outdir, c):
    cl_out = os.path.join(outdir, 'cl_out')
    ttable = os.path.join(outdir, '../truth.txt')
    # this target consumes truth table and clustering output
    source = [ttable, cl_out]
    # this target creates 3 output files
    target = ['{0}.{1}'.format(cl_out, suffix) for suffix in ['joined', 'f1', 'vm']]
    action = 'bin/pbsrun_SCORE.sh $SOURCES.abspath'
    return env.Command(target, source, action)


wrap.add_controls(Environment())
