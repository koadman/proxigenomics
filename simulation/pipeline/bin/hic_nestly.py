from nestly import Nest, stripext
import glob
import os
import os.path

nest = Nest()

commPaths = [os.path.abspath(dn) for dn in glob.glob('references/*')]

# Variation
nest.add('community', commPaths, label_func=os.path.basename)
nest.add('wgs_xfold', (10, 20, 50))
nest.add('hic_table', ['uniform.table'], label_func=stripext)
nest.add('hic_n_frag', (100000, 250000, 500000))

# Constants
nest.add('seed', [2136841], create_dir=False)
nest.add('refseq', ['genomes.fasta'], create_dir=False)
nest.add('base_dir', [os.getcwd()], create_dir=False)
nest.add('wgs_read_length', [150], create_dir=False)
nest.add('wgs_insert_length', [450], create_dir=False)
nest.add('wgs_insert_sd', [100], create_dir=False)
nest.add('wgs_base', ['wgs'],  create_dir=False)
nest.add('hic_base', ['hic'], create_dir=False)
nest.add('hic_inter_prob', [0.9], create_dir=False)
nest.add('hic_read_length', [150], create_dir=False)
#nest.add('hic_n_frag', [250000], create_dir=False)
nest.add('scf2ref', ['scf2ref'], create_dir=False)
nest.add('wgs2scf', ['wgs2scf'], create_dir=False)
nest.add('hic2scf', ['hic2scf'], create_dir=False)

nest.build('runs1')