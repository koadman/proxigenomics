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

wrap.add('lognorm_rel_abundance_mu', config['reference']['lognorm_rel_abundance_mu'])
wrap.add('lognorm_rel_abundance_sigma', config['reference']['lognorm_rel_abundance_sigma'])
wrap.add('num_samples', config['reference']['samples'])

wrap.add('wgs_xfold', config['wgs_xfold'])


@wrap.add_target('make_wgs2ctg')
@name_targets
def make_wgs2ctg(outdir, c):
    com = c['community']
    mu = c['lognorm_rel_abundance_mu']
    sigma = c['lognorm_rel_abundance_sigma']

    for i in range(0,c['num_samples']):
        # TODO find a better way to obtain the path to WGS reads
        query = appconfig.get_wgs_reads_by_sample(
                    os.path.join(os.path.abspath(config['wgs_folder']),
                    com, mu, sigma, str(c['wgs_xfold'])), str(c['num_samples']),
                    config)

        subject = os.path.join(os.path.abspath(config['wgs_folder']),
                    com, mu, sigma, str(c['num_samples']), str(c['wgs_xfold'])), config['wgs_asmdir'],
                    '{0[wgs_base]}.contigs.fasta'.format(config))

        target = os.path.join(outdir, config['wgs2ctg'])
        source = [subject] + query

        action = exec_env.resolve_action({
            'pbs': 'bin/pbsrun_BWA.sh $SOURCES.abspath $TARGET.abspath',
            'sge': 'bin/sgerun_BWA.sh $SOURCES.abspath $TARGET.abspath'
        })
        env.Command(target, source, action)

    return 'output', 0

@wrap.add_target('make_deconvolve')
@name_targets
def make_deconvolve(outdir, c):
    com = c['community']
    mu = c['lognorm_rel_abundance_mu']
    sigma = c['lognorm_rel_abundance_sigma']

    subject = os.path.join(os.path.abspath(config['wgs_folder']),
                com, mu, sigma, str(c['num_samples']), str(c['wgs_xfold'])), config['wgs_asmdir'],
                '{0[wgs_base]}.contigs.fasta'.format(config))


    query_sequences = ""
    for i in range(0,c['num_samples']):
        # TODO find a better way to obtain the path to WGS reads
        query = appconfig.get_wgs_reads_by_sample(
                    os.path.join(os.path.abspath(config['wgs_folder']),
                    com, mu, sigma, str(c['wgs_xfold'])), str(c['num_samples']),
                    config)
        query_sequences = query_sequences + " " + query

    target = os.path.join(outdir, config['wgs2ctg'])
    source = [subject] + query_sequences

    action = exec_env.resolve_action({
        'pbs': 'bin/pbsrun_DECONVOLVE.sh $SOURCES.abspath $TARGET.abspath',
        'sge': 'bin/sgerun_DECONVOLVE.sh $SOURCES.abspath $TARGET.abspath'
    })
    env.Command(target, source, action)

    return 'output', 0


wrap.add_controls(Environment())
