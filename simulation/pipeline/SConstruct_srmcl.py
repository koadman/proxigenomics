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
                                    config['cluster']['algorithms']['srmcl']['folder']))
env = Environment(ENV=os.environ)

# Used for resolving what type of execution environment will be used.
exec_env = appconfig.ExecutionEnvironment(ARGUMENTS, supported_env=['pbs', 'sge', 'local'])

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
        'pbs': 'bin/pbsrun_GRAPH.sh -s $SOURCE.abspath $TARGETS.abspath',
        'sge': 'bin/sgerun_GRAPH.sh -s $SOURCE.abspath $TARGETS.abspath',
        'local': 'bin/bamToEdges.py -s $SOURCE.abspath $TARGETS.abspath'
    })

    return 'edges', 'nodes', env.Command(targets, source, action)


@wrap.add_target('make_cluster_input')
@name_targets
def make_cluster_input(outdir, c):

    sources = [str(c['make_graph']['edges']), str(c['make_graph']['nodes'])]
    base_out = appconfig.prepend_paths(outdir, config['cluster']['input'])[0]
    targets = [base_out, base_out + '.nodemap']

    action = exec_env.resolve_action({
        'pbs': 'bin/pbsrun_MKMETIS.sh {0[ctg_minlen]} $SOURCES.abspath $TARGETS.abspath'.format(config),
        'sge': 'bin/sgerun_MKMETIS.sh {0[ctg_minlen]} $SOURCES.abspath $TARGETS.abspath'.format(config),
        'local': 'bin/edgeToMetis.py --fmt metis -m {0[ctg_minlen]} $SOURCES.abspath $TARGETS.abspath'.format(config)
    })

    return 'output', 'nodemap', env.Command(targets, sources, action)

params = config['cluster']['algorithms']['srmcl']
wrap.add('inflation', numpy.linspace(params['infl']['min'], params['infl']['max'], params['infl']['steps']))
#wrap.add('balance', numpy.linspace(params['bal']['min'], params['bal']['max'], params['bal']['steps']))
# These are added for future sweep possibility. Defaults for now
wrap.add('balance', [0.5])
wrap.add('penalty', [1.25])
wrap.add('redundancy', [0.6])
wrap.add('quality', [0])

@wrap.add_target('do_cluster')
@name_targets
def do_cluster(outdir, c):
    # TODO run over both weighted/unweighted?
    source = c['make_cluster_input']['output']
    target = appconfig.prepend_paths(outdir, config['cluster']['output'] + '.metis')

    action = exec_env.resolve_action({
        'pbs': 'bin/pbsrun_SRMCL.sh -b {0[balance]} -i {0[inflation]} $SOURCE.abspath $TARGET.abspath'.format(c),
        'sge': 'bin/sgerun_SRMCL.sh -b {0[balance]} -i {0[inflation]} $SOURCE.abspath $TARGET.abspath'.format(c),
        'local': 'bin/srmcl -b {0[balance]} -i {0[inflation]} -o $TARGET.abspath $SOURCE.abspath'.format(c)
    })

    return 'output', env.Command(target, source, action)


@wrap.add_target('do_convert')
@name_targets
def do_cluster(outdir, c):
    # TODO run over both weighted/unweighted?
    sources = [c['make_cluster_input']['nodemap'], c['do_cluster']['output']]
    target = appconfig.prepend_paths(outdir, config['cluster']['output'])

    action = exec_env.resolve_action({
        'pbs': 'bin/metisClToMCL.py $SOURCES.abspath $TARGET.abspath'.format(c),
        'sge': 'bin/metisClToMCL.py $SOURCES.abspath $TARGET.abspath'.format(c),
        'local': 'bin/metisClToMCL.py $SOURCES.abspath $TARGET.abspath'.format(c)
    })

    return 'output', env.Command(target, sources, action)


@wrap.add_target('do_score')
def do_score(outdir, c):
    cl_out = c['do_convert']['output']

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
        'sge': 'bin/sgerun_SCORE.sh $SOURCES.abspath',
        'local': 'bin/all_scores.sh $SOURCES.abspath'
    })

    return env.Command(target, source, action)


wrap.add_controls(Environment())
