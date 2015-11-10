from nestly import Nest, stripext
from nestly.scons import SConsWrap, name_targets
import os
import os.path
import appconfig

config = appconfig.read('config.yaml')

nest = Nest()
wrap = SConsWrap(nest, config['wgs_folder'])
env = Environment(ENV=os.environ)

# Constants
wrap.add('seed', [config['seed']], create_dir=False)
wrap.add('refseq', [config['community']['seq']], create_dir=False)
wrap.add('wgs_read_length', [config['wgs_read_length']], create_dir=False)
wrap.add('wgs_insert_length', [config['wgs_insert_length']], create_dir=False)
wrap.add('wgs_insert_sd', [config['wgs_insert_sd']], create_dir=False)
wrap.add('wgs_base', [config['wgs_base']],  create_dir=False)

# Variation
genomes = appconfig.find_files(config['community']['folder'], config['community']['seq'])
commPaths = [os.path.dirname(pn) for pn in genomes]
# For testing - constraint on branch length
#commPaths = [pn for pn in commPaths if float(os.path.basename(pn)) > 0.4 and float(os.path.basename(pn)) < 1]
wrap.add('community', commPaths)

wrap.add('lognorm_rel_abundance_mu', [config['reference']['lognorm_rel_abundance_mu']], create_dir=False)
wrap.add('lognorm_rel_abundance_sigma', [config['reference']['lognorm_rel_abundance_sigma']], create_dir=False)
wrap.add('num_samples', config['reference']['samples'])

wrap.add('wgs_xfold', config['wgs_xfold'])

@wrap.add_target('generate_wgs')
@name_targets
def generate_wgs(outdir, c):
    target = appconfig.get_wgs_reads(outdir, config)
    source = '{1[community][folder]}/{0[community]}/{0[refseq]}'.format(c, config)
    action = 'bin/sgerun_METAART.sh {0[seed]} {0[wgs_insert_length]} {0[wgs_insert_sd]} ' \
             '{0[wgs_read_length]} {0[wgs_xfold]} {0[num_samples]} $SOURCE.abspath {0[wgs_base]} {od} ' \
             '{0[lognorm_rel_abundance_mu]} {0[lognorm_rel_abundance_sigma]} '.format(c, od=outdir)
    for i in range(0,c['num_samples']):
        query = appconfig.get_wgs_reads_by_sample(
             os.path.join(os.path.abspath(config['wgs_folder']),
             c['community'], str(c['num_samples']), str(c['wgs_xfold'])), i,
             config)
        target.append(query)

    return 'abc',env.Command(target, source, action)


#@wrap.add_target('assemble_wgs')
@name_targets
def assemble_wgs(work_dir, c):
    asm_dir = os.path.join(work_dir, config['wgs_asmdir'])
    target = '{od}/{0[wgs_base]}.contigs.fasta'.format(c, od=asm_dir)
    source = [str(c['generate_wgs']['r1']), str(c['generate_wgs']['r2'])]
    action = 'bin/sge_a5submit.sh -mwf -t {0[wgs_base]} $SOURCES.abspath {od}'.format(c, od=asm_dir)

    return 'ctg', env.Command(target, source, action)


@wrap.add_target('index_ctg')
def index_ctg(outdir, c):
    seqFolder = os.path.join(config['community']['folder'], c['community'])
    seq_dat = os.path.join(seqFolder, "ancestral.dat")
    seq_fa = os.path.join(seqFolder, "ancestral.fasta")
    os.system("echo \">ancestor\" > " + seq_fa) 
    os.system("tail -n 1 " + seq_dat + " | wrap >> " + seq_fa) 
    source = seq_fa
    target = source + '.bwt'
    action = 'bin/sgerun_INDEX.sh $SOURCE.abspath'.format(c)

    return env.Command(target, source, action)


wrap.add_controls(Environment())
