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

wrap.add('lognorm_rel_abundance_mu', [config['reference']['lognorm_rel_abundance_mu']], create_dir=False)
wrap.add('lognorm_rel_abundance_sigma', [config['reference']['lognorm_rel_abundance_sigma']], create_dir=False)
wrap.add('num_samples', config['reference']['samples'])

wrap.add('wgs_xfold', config['wgs_xfold'])


@wrap.add_target('make_readmap')
@name_targets
def make_readmap(outdir, c):
    com = c['community']
    mu = c['lognorm_rel_abundance_mu']
    sigma = c['lognorm_rel_abundance_sigma']

    result = ('0','0')
    for i in range(0,c['num_samples']):
        # TODO find a better way to obtain the path to WGS reads
        query = appconfig.get_wgs_reads_by_sample(
                    os.path.join(os.path.abspath(config['wgs_folder']),
                    com, str(c['num_samples']), str(c['wgs_xfold'])), i,
                    config)

        subject = os.path.join(os.path.abspath(config['wgs_folder']),
                    com, str(c['num_samples']), str(c['wgs_xfold']), config['wgs_asmdir'],
                    '{0[wgs_base]}.contigs.fasta'.format(config))

        target = appconfig.get_bam_by_sample(outdir, i, config)
        source = [subject] + query

        action = exec_env.resolve_action({
            'pbs': 'bin/pbsrun_BWA.sh $SOURCES.abspath $TARGET.abspath',
            'sge': 'bin/sgerun_BWA.sh $SOURCES.abspath $TARGET.abspath'
        })

        result = env.Command(target, source, action)
    return 'bams',result

@wrap.add_target('make_deconvolve')
@name_targets
def make_deconvolve(outdir, c):
    com = c['community']

    subject = os.path.join(os.path.abspath(config['wgs_folder']),
                com, str(c['num_samples']), str(c['wgs_xfold']), config['wgs_asmdir'],
                '{0[wgs_base]}.contigs.fasta'.format(config))


    bam_files = []
    for i in range(0,c['num_samples']):
        # TODO find a better way to obtain the path to WGS reads
        bam = appconfig.get_bam_by_sample(
                    os.path.join(os.path.abspath(config['map_folder']),
                    com, str(c['num_samples']),str(c['wgs_xfold'])), i,
                    config)
        bam_files = bam_files + [bam]

    target = os.path.join(outdir, "strains.tre")
    source = [subject] + [bam_files]

    action = exec_env.resolve_action({
        'pbs': 'bin/pbsrun_DECONVOLVE.sh ' + outdir + ' 4 $SOURCES.abspath $TARGET.abspath',
        'sge': 'bin/sgerun_DECONVOLVE.sh ' + outdir + ' 4 $SOURCES.abspath $TARGET.abspath'
    })
    return 'tree',env.Command(target, source, action)



wrap.add_controls(Environment())
