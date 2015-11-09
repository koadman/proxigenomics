#!/usr/bin/env perl

# This script automates the process of generating a set of evolved
# sequences using specific evolution parameters for seq-gen and
# sgEvolver, aligning those sequences, and scoring the alignments


use strict;
use POSIX;
require simujobparams;

if( @ARGV != 2 ){
	die "Usage: simujobrun.pl <raw input sequence> <random seed>";
}

$simujobparams::ancestral_donor = getcwd()."/ancestral_donor.raw";
$simujobparams::seqgen_random = $ARGV[1];
$simujobparams::seqgen_donor_random = $ARGV[1]+1;
$simujobparams::sgEvolver_random = $ARGV[1]+2;

# read in the ancestral sequence FastA and make it raw
my $ancestor_length = 0;
open(INFASTA, $ARGV[0]) || die "Unable to read input sequence $ARGV[0]\n";
open(OUTRAW, ">$simujobparams::ancestral_donor");
while(my $line = <INFASTA>){
	chomp $line;	
	next if $line =~ /^>/;
	$line =~ s/[^ACGT]/A/g;
	print OUTRAW $line;
	$ancestor_length += length($line);
}
close OUTRAW;
close INFASTA;

# fix the sequence length to the complete ancestral genome
if ($simujobparams::seq_length > $ancestor_length/1.5) {
    die "Requested evolved sequence length too long in comparison to input sequence"
}
# ancestor from the beginning of input sequence
$simujobparams::ancestral_start = 0;
$simujobparams::ancestral_length = $simujobparams::seq_length;
# donor uses remainder of input sequence
$simujobparams::donor_start = $simujobparams::seq_length;
$simujobparams::donor_length = $ancestor_length - $simujobparams::seq_length;

# Set the location of the evolution and scoring tools here:
# this is the directory where seq-gen, sgEvolver, and scoreAlignment reside
my $xmfa_alignment = "";

# an array of files that need to be deleted when cleaning up
my @delete_files = ();


# ROADMAP:
# 1) simulate evolution according to parameters
# 2) align genomes using the selected aligner
# 3) score alignments


		# write any command lines executed to a file
open( COMLINES, ">command_lines.txt"  );


		# generate an ancestral base
open( ANCESTRAL, ">$simujobparams::ancestral_seq_name" );
print ANCESTRAL " 1 $simujobparams::seq_length\n";
print ANCESTRAL "Ancestral\n";
close ANCESTRAL;
	
my $rval;

		# extract the specified ancestral sequence
my $extract_cl = "dd if=$simujobparams::ancestral_donor bs=1 skip=$simujobparams::ancestral_start count=$simujobparams::seq_length";
executeCommand( $extract_cl, ">$simujobparams::ancestral_seq_name", "ancestral_dd.err" );

		# remove any stale output file
my $rm_cl = "rm -f $simujobparams::seqgen_out_name";
executeCommand( $rm_cl, "rm.out", "rm.err" );

		# generate a data set
my $seqgen_cl = $simujobparams::tools_dir."seq-gen -mHKY -a $simujobparams::gamma_shape -i 0 -z $simujobparams::seqgen_random -f ";
$seqgen_cl .= "$simujobparams::nt_a_freq $simujobparams::nt_c_freq $simujobparams::nt_g_freq $simujobparams::nt_t_freq ";
$seqgen_cl .= "-t 2 -k 1 -wa -on -s $simujobparams::nt_sub_scale $simujobparams::tree_filename ";
$seqgen_cl .= "< $simujobparams::ancestral_seq_name ";
executeCommand( $seqgen_cl, ">$simujobparams::seqgen_out_name", "seqgen_ancestral.err" );

		# generate a donor base
open( DONOR, ">$simujobparams::donor_seq_name" );
print DONOR " 1 $simujobparams::donor_length\n";
print DONOR "Donor\n";
close DONOR;


		# extract the specified donor sequence
$extract_cl = "dd if=$simujobparams::ancestral_donor bs=1 skip=$simujobparams::donor_start count=$simujobparams::donor_length";
executeCommand( $extract_cl, ">$simujobparams::donor_seq_name", "donor_dd.err" );

		# generate a donor data set
$seqgen_cl = $simujobparams::tools_dir."seq-gen -mHKY -a $simujobparams::gamma_shape -i 0 -z $simujobparams::seqgen_donor_random  -f ";
$seqgen_cl .= "$simujobparams::nt_a_freq $simujobparams::nt_c_freq $simujobparams::nt_g_freq $simujobparams::nt_t_freq ";
$seqgen_cl .= "-t 2 -k 1 -wa -on -s $simujobparams::nt_sub_scale $simujobparams::tree_filename ";
$seqgen_cl .= "< $simujobparams::donor_seq_name ";
executeCommand( $seqgen_cl, ">$simujobparams::seqgen_out_name", "seqgen_donor.err" );


		# continue evolution with sgEvolver
my $sgevolver_cl = $simujobparams::tools_dir."sgEvolver --indel-freq=".$simujobparams::indel_rate." --small-ht-freq=$simujobparams::small_ht_rate --small-ht-size=$simujobparams::small_ht_size".
" --large-ht-freq=$simujobparams::large_ht_rate --inversion-freq=$simujobparams::inv_rate --large-ht-min=$simujobparams::large_ht_min --large-ht-max=$simujobparams::large_ht_max".
" --random-seed=$simujobparams::sgEvolver_random".
" --inversion-size=$simujobparams::inv_size $simujobparams::tree_filename $simujobparams::seqgen_out_name $simujobparams::evolved_seqs_name $simujobparams::evolved_seqs_fname";
$rval = executeCommand( $sgevolver_cl, "sgEvolver.out", "sgEvolver.err" );

		# die if sgEvolver failed
die "Failure in sgEvolver" if( $rval != 0 );

die "Failure in sgEvolver" unless -e "$simujobparams::evolved_seqs_name";
die "Failure in sgEvolver" unless -e "$simujobparams::evolved_seqs_fname";


		# delete seq-gen related files
push( @delete_files, $simujobparams::seqgen_out_name );
#push( @delete_files, $simujobparams::ancestral_seq_name );
push( @delete_files, $simujobparams::donor_seq_name );
		# delete evolved sequence data
push( @delete_files, $simujobparams::evolved_seqs_name );
deleteFiles( @delete_files );
exit(0);

sub executeCommand {
  my $command = shift;
  my $stdout_file = shift;
  my $stderr_file = shift;
  $command .= " >$stdout_file 2>$stderr_file";
  print "Executing $command\n";
  print COMLINES "$command\n";
  my $rval = system($command);
  `echo "Exited with code $rval" >> $stderr_file` if $rval != 0;
  return $rval;
}

sub deleteFiles {
  if( $debug != 1 ){
    foreach(@_){
      `rm -f $_`;
      print COMLINES "rm -f $_\n";
    }
  }
}


