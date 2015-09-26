from nestly import Nest
from nestly.scons import SConsWrap, name_targets
import os
import os.path
import numpy
import appconfig


config = appconfig.read('config.yaml')

# Base folder
nest = Nest()
wrap = SConsWrap(nest, os.path.join(config['cluster']['folder'],
                                    config['cluster']['algorithms']['mcl']['folder']))
env = Environment(ENV=os.environ)

# Used for resolving what type of execution environment will be used.
exec_env = appconfig.ExecutionEnvironment(ARGUMENTS, supported_env=['pbs', 'local'])

# Variation

# don't include root as we don't want it embedded in this nest hierarchy
hic_paths = appconfig.get_precedents(config['map_folder'], config['hic2ctg'], prepend_root=False)
wrap.add('hic_path', hic_paths)


#
# TODO, this needs to be placed in mapping
#
@wrap.add_target('make_graph')
@name_targets
def make_graph(outdir, c):
    # add the root back in because we need to refer to the file
    ref_path = os.path.join(config['map_folder'], c['hic_path'])
    hic_bam = str(os.path.join(ref_path, config['hic2ctg']))

    source = hic_bam
    targets = appconfig.prepend_paths(outdir, ['edges.csv', 'nodes.csv'])
    action = exec_env.resolve_action({
        'pbs': 'bin/pbsrun_GRAPH.sh $SOURCE.abspath $TARGETS.abspath',
        'local': 'bin/bamToEdges.py $SOURCE.abspath $TARGETS.abspath'
    })

    return 'edges', 'nodes', env.Command(targets, source, action)


@wrap.add_target('make_cluster_input')
@name_targets
def make_cluster_input(outdir, c):

    source = [str(c['make_graph']['edges']), str(c['make_graph']['nodes'])]
    target = appconfig.prepend_paths(outdir, config['cluster']['input'])

    action = exec_env.resolve_action({
        'pbs': 'bin/pbsrun_MKMCL.sh {0[ctg_minlen]} $SOURCES.abspath $TARGET.abspath'.format(config),
        'local': 'bin/makeMCLinput.py {0[ctg_minlen]} $SOURCES.abspath $TARGET.abspath'.format(config)
    })

    return 'output', env.Command(target, source, action)

mcl_infl = config['cluster']['algorithms']['mcl']['infl']
wrap.add('mcl_inflation', numpy.linspace(mcl_infl['min'], mcl_infl['max'], mcl_infl['steps']))

@wrap.add_target('do_cluster')
@name_targets
def do_cluster(outdir, c):
    # TODO run over both weighted/unweighted?
    source = c['make_cluster_input']['output']
    target = appconfig.prepend_paths(outdir, config['cluster']['output'])

    action = exec_env.resolve_action({
        'pbs': 'bin/pbsrun_MCL.sh {0[mcl_inflation]} $SOURCE.abspath $TARGET.abspath'.format(c),
        'local': 'bin/mcl $SOURCE.abspath --abc -I {0[mcl_inflation]} -o $TARGET.abspath'.format(c)
    })

    return 'output', env.Command(target, source, action)


@wrap.add_target('do_score')
def do_score(outdir, c):
    cl_out = c['do_cluster']['output']

    ref_path = os.path.join(config['map_folder'], c['hic_path'])
    ttable = appconfig.search_up(ref_path, config['truth_table'])
    if ttable is None:
        raise RuntimeError('Could not find an accompanying truth table for associated run {0}'.format(c['hic_path']))

    # this target consumes truth table and clustering output
    source = [ttable, cl_out]
    # this target creates 3 output files
    target = ['{0}.{1}'.format(cl_out, suffix) for suffix in ['f1', 'vm', 'bc']]

    action = exec_env.resolve_action({
        'pbs': 'bin/pbsrun_SCORE.sh $SOURCES.abspath',
        'local': 'bin/all_scores.sh $SOURCES.abspath'
    })

    return env.Command(target, source, action)


wrap.add_controls(Environment())
