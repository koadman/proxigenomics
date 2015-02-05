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
    return os.path.join(os.path.abspath('wgs_data'), com, tbl,
                        str(c['wgs_xfold']), 'wgs', 'wgs.contigs.fasta')


def prepend_path(path, fnames):
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
wrap.add_aggregate('mcl_files', list)
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

    action = 'bin/pbsrun_MAP.sh $SOURCES.abspath $TARGET.abspath'.format(od=outdir)

    c['sam_files'].append(target)

    return env.Command(target, source, action)


@wrap.add_target('make_ctg2ref')
def map_ctg2ref(outdir, c):
    com = os.path.basename(c['community'])
    tbl = stripext(c['hic_table'])

    query = os.path.join(
                os.path.abspath(config['wgs_folder']),
                com, tbl, str(c['wgs_xfold']), config['wgs_asmdir'],
                '{0[wgs_base]}.contigs.fasta'.format(config))

    subject = os.path.join(c['community'], config['community']['seq'])

    target = os.path.join(outdir, config['ctg2ref'])

    source = [subject, query]

    action = 'bin/pbsrun_MAP.sh $SOURCES.abspath $TARGET.abspath'.format(od=outdir)

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
    action = 'bin/pbsrun_MAP2.sh $SOURCES.abspath $TARGET.abspath'.format(od=outdir)

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
    wgs_bam = [os.path.splitext(fn)[0] + ".bam" for fn in c['sam_files'] if fn.endswith('wgs2ctg.sam')]
    source = hic_sam + wgs_bam
    target = prepend_path(outdir, ['edges.csv', 'nodes.csv'])
    c['graph_files'].append(target)
    action = 'cd {od} && /panfs/panspermia/120274/work/hi-c/sim/bin/pbsrun_GRAPH.sh ' \
             '$SOURCES.file $TARGETS.file'.format(od=outdir)
    return env.Command(target, source, action)


wrap.add('cluster_method', ['mcl'])
wrap.add('score_min_length', [1000])


@wrap.add_target('make_truth')
def make_truth(outdir, c):
    seq = get_wgs_fasta(c)
    source = [seq] + pick_sam('ctg2ref.sam', c)
    target = os.path.join(outdir, "truth.txt")
    c['truth'].append(target)
    action = 'cd {od} && /panfs/panspermia/120274/work/hi-c/sim/bin/parseSamCigar.py {0[score_min_length]} $SOURCES.abspath $TARGET.file'.format(
        c, od=outdir)
    return env.Command(target, source, action)


@wrap.add_target('make_mclinput')
def make_mclinput(outdir, c):
    source = c['graph_files']
    target = prepend_path(outdir, ['mclIn.weighted', 'mclIn.unweighted'])
    action = 'cd {od} && /panfs/panspermia/120274/work/hi-c/sim/bin/makeMCLinput.py {0[score_min_length]} $SOURCES.abspath'.format(
        c, od=outdir)
    c['mcl_files'].append(target)
    return env.Command(target, source, action)


wrap.add('mcl_inflation', numpy.linspace(1.1, 2.0, 19))


@wrap.add_target('do_mcl')
def do_mcl(outdir, c):
    # figure out how to run over both weighted/unweighted
    source = c['mcl_files'][0][0]
    target = prepend_path(outdir, ['mclout'])
    action = 'cd {od} && /panfs/panspermia/120274/work/hi-c/sim/bin/pbsrun_MCL.sh {0[mcl_inflation]} $SOURCE.abspath $TARGET.file'.format(
        c, od=outdir)
    return env.Command(target, source, action)


@wrap.add_target('do_mcljointruth')
def do_mcl(outdir, c):
    source = c['truth'] + prepend_path(outdir, ['mclout'])
    target = prepend_path(outdir, ['mcl.joined.truth'])
    action = 'cd {od} && /panfs/panspermia/120274/work/hi-c/sim/bin/mclJoinWithTruth.py $SOURCES.abspath $TARGET.file'.format(
        c, od=outdir)
    return env.Command(target, source, action)


@wrap.add_target('do_f1')
def do_score(outdir, c):
    source = os.path.join(outdir, 'mcl.joined.truth')
    target = os.path.join(outdir, 'f1')
    action = 'cd {od} && /panfs/panspermia/120274/work/hi-c/sim/bin/f1score.py $SOURCE.file $TARGET.file'.format(c,
                                                                                                                 od=outdir)
    return env.Command(target, source, action)


@wrap.add_target('do_vmeasure')
def do_vmeasure(outdir, c):
    source = os.path.join(outdir, 'mcl.joined.truth')
    target = os.path.join(outdir, 'vmeasure')
    action = 'cd {od} && /panfs/panspermia/120274/work/hi-c/sim/bin/vmeasure.py $SOURCE.file $TARGET.file'.format(c,
                                                                                                                  od=outdir)
    return env.Command(target, source, action)


wrap.add_controls(Environment())
