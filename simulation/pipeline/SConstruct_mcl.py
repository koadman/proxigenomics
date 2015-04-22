from nestly import Nest, stripext
from nestly.scons import SConsWrap, name_targets
import os
import os.path
import numpy
import appconfig
import types

#
# Helper functions
#

def prepend_paths(path, fnames):
    if isinstance(fnames, types.StringTypes):
        fnames = [fnames]
    return [os.path.join(path, fn) for fn in fnames]


config = appconfig.read('config.yaml')

# Base folder
nest = Nest()
wrap = SConsWrap(nest, config['cluster']['folder'])
env = Environment(ENV=os.environ)

# Variation
hic_bams = appconfig.find_files(config['map_folder'], config['hic2ctg'])
hic_paths = [os.path.join(config['map_folder'], os.path.dirname(pn)) for pn in hic_bams]
wrap.add('hic_path', hic_paths)

@wrap.add_target('make_input')
@name_targets
def make_graph(outdir, c):
    hic_bam = str(os.path.join(c['hic_path'], config['hic2ctg']))
    wgs_bam = appconfig.search_up(config['map_folder'], c['hic_path'], config['wgs2ctg'])
    if wgs_bam is None:
        raise RuntimeError('Could not find an accompanying wgs bam for hic bam {0}'.format(hic_bam))

    sources = [hic_bam, wgs_bam]
    target = prepend_paths(outdir, ['edges.csv', 'nodes.csv'])
    action = 'bin/pbsrun_GRAPH.sh $SOURCES.abspath $TARGETS.abspath'
    return 'edges', 'nodes', env.Command(target, sources, action)


'''
@wrap.add_target('make_cluster_input')
@name_targets
def make_cluster_input(outdir, c):

    source = [str(c['make_graph']['edges']), str(c['make_graph']['nodes'])]
    target = prepend_paths(outdir, config['cluster']['input'])

    if c['cluster_method'] == 'mcl':
        action = 'bin/pbsrun_MKMCL.sh {1[ctg_minlen]} $SOURCES.abspath $TARGET.abspath'.format(c, config)

    elif c['cluster_method'] == 'oclustr':
        action = 'bin/edgeToMetis.py -m {1[ctg_minlen]} -f graphml $SOURCES.abspath $TARGET.abspath'.format(c, config)

    return 'output', env.Command(target, source, action)

mcl_param = config['cluster']['mcl_infl']
wrap.add('mcl_inflation', numpy.linspace(mcl_param['min'], mcl_param['max'], mcl_param['steps']))

@wrap.add_target('do_mcl')
@name_targets
def do_mcl(outdir, c):
    # TODO run over both weighted/unweighted?
    source = c['make_cluster_input']['output']
    target = prepend_paths(outdir, config['cluster']['output'])

    action = 'bin/pbsrun_MCL.sh {0[mcl_inflation]} $SOURCE.abspath $TARGET.abspath'.format(c)

    return 'output', env.Command(target, source, action)


@wrap.add_target('do_score')
def do_score(outdir, c):
    #cl_out = os.path.join(outdir, config['cluster']['output'])
    cl_out = c['do_mcl']['output']
    ttable = c['make_truth']['output']
    # this target consumes truth table and clustering output
    source = [ttable, cl_out]
    # this target creates 3 output files
    target = ['{0}.{1}'.format(cl_out, suffix) for suffix in ['f1', 'vm', 'bc']]

    action = 'bin/pbsrun_SCORE.sh $SOURCES.abspath'

    return env.Command(target, source, action)
'''

wrap.add_controls(Environment())
