from nestly import Nest
from nestly.scons import SConsWrap, name_targets
import os
import os.path
import appconfig


config = appconfig.read('config.yaml')

# Base folder
nest = Nest()
wrap = SConsWrap(nest, os.path.join(config['cluster']['folder'],
                                    config['cluster']['algorithms']['louvain']['folder']))
env = Environment(ENV=os.environ)

# Used for resolving what type of execution environment will be used.
exec_env = appconfig.ExecutionEnvironment(ARGUMENTS)

# don't include root as we don't want it embedded in this nest hierarchy
hic_paths = appconfig.get_precedents(config['map_folder'], config['hic2ctg'], prepend_root=False)
wrap.add('hic_path', hic_paths)

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
        'pbs': 'bin/edgeToMetis.py -m {0[ctg_minlen]} -f graphml $SOURCES.abspath $TARGET.abspath'.format(config),
        'local': 'bin/edgeToMetis.py -m {0[ctg_minlen]} -f graphml $SOURCES.abspath $TARGET.abspath'.format(config)
    })

    return 'output', env.Command(target, source, action)


wrap.add('otype', config['cluster']['algorithms']['louvain']['otype'])

@wrap.add_target('do_cluster')
@name_targets
def do_cluster(outdir, c):
    # TODO run over both weighted/unweighted?
    source = c['make_cluster_input']['output']
    target = appconfig.prepend_paths(outdir, config['cluster']['output'])

    action = exec_env.resolve_action({
        'pbs': 'bin/pbsrun_LOUVAIN.sh {0[otype]} $SOURCE.abspath $TARGET.abspath'.format(c),
        'local': 'bin/louvain.py --otype {0[otype]} $SOURCE.abspath $TARGET.abspath'.format(c)
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

    action =  exec_env.resolve_action({
        'pbs': 'bin/pbsrun_SCORE.sh $SOURCES.abspath',
        'local': 'bin/all_scores.sh $SOURCES.abspath'
    })

    return env.Command(target, source, action)


wrap.add_controls(Environment())
