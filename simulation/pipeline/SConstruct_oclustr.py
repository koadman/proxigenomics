from nestly import Nest
from nestly.scons import SConsWrap, name_targets
import os
import os.path
import appconfig


config = appconfig.read('config.yaml')

# Base folder
nest = Nest()
wrap = SConsWrap(nest, os.path.join(config['cluster']['folder'],
                                    config['cluster']['algorithms']['oclustr']['folder']))
env = Environment(ENV=os.environ)

# control whether local or distributed execution of targets
execution_type = ARGUMENTS.get('exec_type', 'pbs')
if not (execution_type == 'pbs' or execution_type == 'local'):
    raise RuntimeError('Unknown execution type specified. Select either "pbs" or "local" [default: pbs]')

def resolve_action(action_dict):
    """
    Simple method to resolve the action associated with a particular
    execution environment.
    :param action_dict: the set of actions to choose from.
    :return: the action associated with the value of global execution_type
    """
    return action_dict[execution_type]

# don't include root as we don't want it embedded in this nest hierarchy
hic_paths = appconfig.get_precedents(config['map_folder'], config['hic2ctg'], prepend_root=False)
wrap.add('hic_path', hic_paths)

@wrap.add_target('make_graph')
@name_targets
def make_graph(outdir, c):
    # add the root back in because we need to refer to the file
    ref_path = os.path.join(config['map_folder'], c['hic_path'])
    hic_bam = str(os.path.join(ref_path, config['hic2ctg']))
    wgs_bam = appconfig.search_up(ref_path, config['wgs2ctg'])
    if wgs_bam is None:
        raise RuntimeError('Could not find an accompanying wgs bam for ref path {0}'.format(ref_path))

    sources = [hic_bam, wgs_bam]
    target = appconfig.prepend_paths(outdir, ['edges.csv', 'nodes.csv'])

    action = resolve_action({
        'pbs': 'bin/pbsrun_GRAPH.sh $SOURCES.abspath $TARGETS.abspath',
        'local': 'bin/bamToEdges.py --wgs $SOURCES.abspath $TARGETS.abspath'
    })

    return 'edges', 'nodes', env.Command(target, sources, action)


@wrap.add_target('make_cluster_input')
@name_targets
def make_cluster_input(outdir, c):

    source = [str(c['make_graph']['edges']), str(c['make_graph']['nodes'])]
    target = appconfig.prepend_paths(outdir, config['cluster']['input'])

    action = resolve_action({
        'pbs': 'bin/edgeToMetis.py -m {0[ctg_minlen]} -f graphml $SOURCES.abspath $TARGET.abspath'.format(config),
        'local': 'bin/edgeToMetis.py -m {0[ctg_minlen]} -f graphml $SOURCES.abspath $TARGET.abspath'.format(config)
    })

    return 'output', env.Command(target, source, action)

wrap.add('isolates', ['isolates', 'no-isolates'])

@wrap.add_target('do_cluster')
@name_targets
def do_mcl(outdir, c):
    # TODO run over both weighted/unweighted?
    source = c['make_cluster_input']['output']
    target = appconfig.prepend_paths(outdir, config['cluster']['output'])

    if c['isolates'] == 'isolates':
        action = resolve_action({
            'pbs': 'bin/pbsrun_OCLUSTR.sh -i $SOURCE.abspath $TARGET.abspath',
            'local': 'bin/oclustr.py $SOURCE.abspath $TARGET.abspath'
        })
    else:
        action = resolve_action({
            'pbs': 'bin/pbsrun_OCLUSTR.sh $SOURCE.abspath $TARGET.abspath',
            'local': 'bin/oclustr.py --no-isolates $SOURCE.abspath $TARGET.abspath'
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

    action = resolve_action({
        'pbs:': 'bin/pbsrun_SCORE.sh $SOURCES.abspath',
        'local': 'bin/all_scores.sh $SOURCES.abspath'
    })

    return env.Command(target, source, action)


wrap.add_controls(Environment())
