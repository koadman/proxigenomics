from nestly import Nest, stripext
import glob
import os
import os.path

nest = Nest()

commPaths = [os.path.abspath(dn) for dn in glob.glob('references/*')]

# Variation
nest.add('community', commPaths, label_func=os.path.basename)
nest.add('hic_table', ['uniform.table'], label_func=stripext)
nest.add('wgs_xfold', [10, 20, 50])

# Constants
nest.add('seed', [2136841], create_dir=False)
nest.add('refseq', ['genomes.fasta'], create_dir=False)
nest.add('base_dir', [os.getcwd()], create_dir=False)
nest.add('wgs_read_length', [150], create_dir=False)
nest.add('wgs_insert_length', [450], create_dir=False)
nest.add('wgs_insert_sd', [100], create_dir=False)
nest.add('wgs_base', ['wgs'],  create_dir=False)

nest.build('test')
