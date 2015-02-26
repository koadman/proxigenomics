#!/usr/bin/env python
from Bio import Phylo
import argparse
import os.path

parser = argparse.ArgumentParser(description='Prepare the parameter file for running sgEvolver')
parser.add_argument('--tree', required=True, nargs=1, help='Phylogenetic tree')
parser.add_argument('--seq', required=True, nargs=1, help='Raw ancestral sequence')
parser.add_argument('--seq-len', required=True, type=int, nargs=1, help='Length of generated sequences (bp)')
parser.add_argument('--tree-scale', default=1, type=float, nargs=1, help="Tree scale factor")
parser.add_argument('--fmt', choices=['newick', 'nexus', 'nexml', 'phyloxml', 'cdao'],
                    default='newick', nargs=1, help='Tree format')
parser.add_argument('-o', dest='output', required=True, nargs=1, help='Output path')
args = parser.parse_args()

if not os.path.exists(args.output[0]) or not os.path.isdir(args.output[0]):
    raise IOError('Output path {0} is not a directory or does not exist'.format(args.output[0]))

out_file = os.path.join(args.output[0], 'simujobparams.pm')

tree = Phylo.read(args.tree[0], args.fmt)
tips = tree.get_terminals()
seq_names = sorted([t.name for t in tips])
seq_count = len(seq_names)

print 'Found {0} tip tables: {1}\n'.format(seq_count, seq_names)
print '\nAscii Tree Representation:'
Phylo.draw_ascii(tree)

with open(out_file, 'w') as h_out:
    h_out.write('package simujobparams;\n\n')
    h_out.write('#\n# user set values\n#\n')
    h_out.write('@seqnames=({0});\n'.format(','.join(seq_names)))
    h_out.write('$seq_count="{0}";\n'.format(seq_count))
    h_out.write('$tree_filename="{0}";\n'.format(args.tree[0]))
    h_out.write('$ancestral_donor="{0}";\n'.format(args.seq[0]))
    h_out.write('$seq_length={0};\n'.format(args.seq_len[0]))
    h_out.write('$tree_scale={0};\n'.format(args.tree_scale[0]))
    h_out.write('#\n# extra constants\n#\n')
    h_out.write('$nt_sub_scale=1;\n' \
                '$tools_dir="";\n' \
                '$mauve_dir="";\n' \
                '$mavid_dir="";\n' \
                '$lagan_dir="";\n' \
                '$tba_dir="";\n' \
                '$gamma_shape=1;\n' \
                '$indel_rate=0.1;\n' \
                '$small_ht_rate=0.05;\n' \
                '$small_ht_size=200;\n' \
                '$large_ht_rate=0.005;\n' \
                '$large_ht_min=10000;\n' \
                '$large_ht_max=60000;\n' \
                '$inv_rate=0.005;\n' \
                '$inv_size=50000;\n' \
                '$nt_a_freq=0.25;\n' \
                '$nt_c_freq=0.25;\n' \
                '$nt_g_freq=0.25;\n' \
                '$nt_t_freq=0.25;\n' \
                '$score_subset="";\n' \
                '$ancestral_seq_name="ancestral.dat";\n' \
                '$donor_seq_name="donor.dat";\n' \
                '$seqgen_out_name="seqgen.dat";\n' \
                '$evolved_seqs_name="evolved.dat";\n' \
                '$evolved_seqs_fname="evolved_seqs.fas";\n' \
                '$ancestral_start=0;\n' \
                '$donor_start=0;\n' \
                '$seqgen_random=1661120051;\n' \
                '$seqgen_donor_random=500482372;\n' \
                '$sgEvolver_random=2127648372;\n')
