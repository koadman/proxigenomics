ART_SOLiD  README (last updated on 11/22/2011) Weichun Huang at whduke@gmail.com                              

DESCRIPTION:

	ART_SOLiD (version 1.0.1) is a simulation program to generate sequence read data of SOLiD sequencing
        reads. ART generates reads according to a SOLiD read error profile. The built-in error profile is an
	empricial error profile summarilized from large SOLiD sequencing data. ART has been using for testing
       	or benchmarking a variety of method or tools for next-generation sequencing data analysis, including 
	read alignment, de novo assembly, detection of SNP, CNV, or other structure variation.
       	
	art_SOLid can generate both single-end and paired-end of SOLiD sequencing platform.  Its outputs include
       	FASTQ read, MAP alignment, and optional SAM alignment files. ALN files can also be converted to UCSC BED
       	files by using aln2bed.pl program included. 

COMPILATION AND INSTALLATION:

	PREREQUISITES: 
			1) GNU g++ 4.0 or above (http://gcc.gnu.org/install) 
			2) GNU gsl library (http://www.gnu.org/s/gsl/) 

	./configure --prefix=$HOME
	make
	make install

EXAMPLES

	In the "examples" subdirectory, the shell script "run_test_examples.sh" gives two examples of using
       	art_SOLiD for read simulation.  To test these two examples, just run the script "run_test_examples.sh"
USAGE

       SINGLE-END READ SIMULATION

          art_SOLiD [options] <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <READ_LEN> <FOLD_COVERAGE>

          Example:
              art_SOLiD seq_reference.fa ./outdir/dat_single_end 32 10

       PAIRED-END READ SIMULATION

          art_SOLiD [ options ] <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <READ_LEN> <FOLD_COVERAGE> <MEAN_FRAG_LEN> <STD_DEV>

          Example:
             art_SOLiD seq_reference.fa ./outdir/dat_paired_end 25 10 500 20

      COMMAND-LINE PARAMETERS

          INPUT_SEQ_FILE     -  the filename of DNA reference sequences in FASTA format
          OUTPUT_FILE_PREFIX -  the prefix or directory of output read data file (*.fq) and read alignment file (*.aln)
	  FOLD_COVERAGE      -  the fold of read coverage over the reference sequences 
	  READ_LEN           -  the length of reads
	  MEAN_FRAG_LEN      -  the average DNA fragment size for paired-end read simulation 
	  STD_DEV            -  the standard deviation of the DNA fragment size for paired-end read simulation 

      OPTIONAL PARAMETERS 
      
         -p specify user's own read profile for simulation. Please the built-in profile as example under the subfolder "SOLiD_profile"
	  
	 
OUTPUT DATA FILES

	*.fq   - FASTQ read data files. For paired-end read simulation, *_R3.fq contains data of the first reads, and *_F3.fq for the second reads.
	*.map  - read alignment files. For paired-end read simulation, *_R3.map has read alignments for the first reads and *_F3.map for the second reads.

	FASTQ file format 

		A FASTQ file contains both sequence bases and quality scores of sequencing reads and is in the following format:  

			@read_id 
			sequence_read 
			+ 
			base_quality_scores 
	
		A base quality score is coded by the ASCII code of a single character, where the quality score is equal to ASCII code of the
       		character minus 33.    

		Example: 
		@1_1_1_F3
	       	T10121011322000311302213102311132
	       	+
	       	163=+48./48<347//,=/84-4)77(''-)
	       	@1_1_2_F3
	       	T01213000012110232021321303011332
	       	+
	       	1>7<-653?01:55./.%3+'61-+40(*)#'

	MAP file format 
		A MAP file has a Header and main Body parts. The header part includes the command used to generate this file and reference
	       	sequence id and length. The header @CM tag for command line, and @SQ for reference sequence.  A header always starts with 
		"##ART" and ends with  "##Header End".

		HEADER EXAMPLE

		##ART_SOLiD     read_length     32
		@CM     ../../bin/MacOS64/art_SOLiD -s ./testSeq.fa ./single_SOLiD_test 32 10 
		@SQ     seq1    7207
	       	@SQ     seq2    3056
		##Header End

		The body part of a MAP file is a tab-delimited text alignment data where each line is the mapping information record of a read. The format of each line
		is below:
		
		ref_seq_id	read_id		ref_map_pos	ref_seq_strand	num_sequencing_errors	err_pos_1  WrongCorrectColor_1	err_pos_2  WrongCorrectColor_2	...	
	
		ref_map_pos - the alignment start position of reference sequence. ref_map_pos is always relative to the strand of reference
       		sequence. That is, ref_map_pos 10 in the plus (+) strand is different from ref_map_pos 10 in the minus (‐) stand.  

		num_sequencing_errors - number of sequencing call errors

		err_pos_n - the zero-indexed read position of the n_th call error
		WrongCorrectColor_n  - the n_th call error representing by two digits number with the 1st number the wrong called color and the 2nd the correct color

		Example: 
		
		seq1	1_1_1_F3	6028	+	4	18	31	22	03	25	32	27	23

	SAM format

       		SAM is a standard format for next-gen sequencing read alignments. The details of the format and examples are available at the links below:	
		1) http://samtools.sourceforge.net/SAM1.pdf
	       	2) http://genome.sph.umich.edu/wiki/SAM	

		Read sequences in a SAM file are in the regular DNA base space. As one change in SOLiD reads in the color space can result multiple changes in base space,
	        so not all mismatches in base space between reads and its reference sequence are of SOLiD sequencing errors.

	BED format

		See the format at UCSC http://genome.ucsc.edu/FAQ/FAQformat.html#format1

		NOTE: both MAP and BED format files use 0-based coordinate system while SAM format uses 1-based coordinate system.
