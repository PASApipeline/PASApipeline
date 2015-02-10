#!/usr/bin/env perl

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");

use strict;
use warnings;
use Getopt::Std;
use Gene_obj;
use Nuc_translator;
use CdbTools;
use Fasta_reader;
use Longest_orf;
use List::Util qw (min max);

use vars qw ($opt_t $opt_m $opt_G $opt_h $opt_v $opt_C $opt_T $opt_R $opt_P $opt_B);

&getopts ('t:m:G:hvCT:R:P:B:');


my $usage =  <<_EOH_;

############################# Options ###############################
#
# Required:
#
# -t transcripts.fasta
#
#  Optional:
#
# -m minimum protein length (default: 100)
#
# -G genetic code (default: universal, options: Euplotes, Tetrahymena, Candida, Acetabularia)
#
# -h print this option menu and quit
# -v verbose
#
# -C complete ORFs only ********
# -B examine both strands (default: only top strand)
# -T top longest ORFs to train Markov Model (hexamer stats) (default: 500)
#
###################### Process Args and Options #####################

_EOH_

    ;


if ($opt_h) {die $usage;}


my $transcripts_file = $opt_t || die $usage;

my $min_prot_length = $opt_m || 100;
my $genetic_code = $opt_G;

my $top_ORFs_train = $opt_T || 500;
my $num_rand_iter = $opt_R || 100;
my $both_strands = $opt_B || 0;

my $COMPLETE_ORFS_ONLY = $opt_C;

$|++;

our $SEE = $opt_v || 0;


my $NO_REPLACE = 0;

if ($genetic_code) {
    &Nuc_translator::use_specified_genetic_code($genetic_code);
}


my $prefix = "longest_orfs";
my $cds_file = "$prefix.cds";
my $gff3_file = "$prefix.gff3";
my $pep_file = "$prefix.pep";


unless ($NO_REPLACE && -s $pep_file && -s $cds_file && -s $gff3_file) {
	
	open (PEP, ">$pep_file") or die $!;
	open (CDS, ">$cds_file") or die $!; 
	open (GFF, ">$gff3_file") or die $!;
	
	
	my $counter = 0;
	
	my $fasta_reader = new Fasta_reader($transcripts_file);
	while (my $seq_obj = $fasta_reader->next()) {
		
		my $acc = $seq_obj->get_accession();
		my $sequence = $seq_obj->get_sequence();
		
		my $longest_orf_finder = new Longest_orf();
		$longest_orf_finder->allow_5prime_partials();
		$longest_orf_finder->allow_3prime_partials();
		
		unless ($both_strands) {
			$longest_orf_finder->forward_strand_only();
		}
		
		my @orf_structs = $longest_orf_finder->capture_all_ORFs($sequence);
		
		@orf_structs = reverse sort {$a->{length}<=>$b->{length}} @orf_structs;
		
		my $orf = shift @orf_structs;
		
		my $start = $orf->{start};
		my $stop = $orf->{stop};
		my $length = int($orf->{length}/3);
		my $orient = $orf->{orient};
		my $protein = $orf->{protein};
		
		if ($length < $min_prot_length) { next; }
		
		my $coords_href = { $start => $stop };
		
		my $gene_obj = new Gene_obj();
		
		$counter++;
		$gene_obj->populate_gene_object($coords_href, $coords_href);
		$gene_obj->{asmbl_id} = $acc;
		
		my $model_id = "m.$counter";
		my $gene_id = "g.$counter";
		
		
		$gene_obj->{TU_feat_name} = $gene_id;
		$gene_obj->{Model_feat_name} = $model_id;
		$gene_obj->{com_name} = "ORF $gene_id $model_id";
		
		my $cds = $gene_obj->create_CDS_sequence(\$sequence);
		
		my $got_start = 0;
		my $got_stop = 0;
		if ($protein =~ /^M/) {
			$got_start = 1;
		} 
		if ($protein =~ /\*$/) {
			$got_stop = 1;
		}
        
		my $prot_type = "";
		if ($got_start && $got_stop) {
			$prot_type = "complete";
		} elsif ($got_start) {
			$prot_type = "3prime_partial";
		} elsif ($got_stop) {
			$prot_type = "5prime_partial";
		} else {
			$prot_type = "internal";
		}
		
		if ($COMPLETE_ORFS_ONLY && $prot_type ne "complete") { next; }
		
		
		print PEP ">$model_id $gene_id type:$prot_type len:$length\n$protein\n";
		
		print CDS ">$model_id $gene_id type:$prot_type len:$length\n$cds\n";
		
		print GFF $gene_obj->to_GFF3_format() . "\n";
		
	}
}

## Train a Markov model based on longest candidate CDS sequences, score all candidates, and select the final set.

# get longest entries
my $top_cds_file = "$cds_file.top_${top_ORFs_train}_longest";
my $cmd = "$FindBin::Bin/../misc_utilities/get_top_longest_fasta_entries.pl $cds_file $top_ORFs_train > $top_cds_file";
&process_cmd($cmd) unless ($NO_REPLACE && -s $top_cds_file);

$cmd = "$FindBin::Bin/../misc_utilities/randomize_sequences.pl $transcripts_file > $transcripts_file.random";
&process_cmd($cmd) unless ($NO_REPLACE && -s "$transcripts_file.random");


# get hexamer scores
$cmd = "$FindBin::Bin/../misc_utilities/seq_n_background_to_logliklihood_vals.pl $top_cds_file $transcripts_file.random > hexamer.scores";
&process_cmd($cmd) unless ($NO_REPLACE && -s "hexamer.scores");

# score all cds entries
$cmd = "$FindBin::Bin/../misc_utilities/score_CDS_liklihood_all_6_frames.pl $cds_file hexamer.scores > $cds_file.scores";
&process_cmd($cmd) unless ($NO_REPLACE && -s "$cds_file.scores");

# get accs for best entries
my $acc_file = "$cds_file.scores.selected";
{
	open (my $ofh, ">$acc_file") or die "Error, cannot write to $acc_file";
	open (my $ifh, "$cds_file.scores") or die "Error, cannot open file $cds_file.scores";
	while (<$ifh>) {
		chomp;
		my ($acc, @scores) = split(/\t/);
		
		my $score_1 = shift @scores;
		my $max_score_other_frame = max(@scores);
		if ($score_1 > 0 && $score_1 > $max_score_other_frame) { 
			print $ofh "$acc\n";
		}
	}
	close $ifh;
	close $ofh;
}

# index the current gff file:
$cmd = "$FindBin::Bin/../misc_utilities/index_gff3_files_by_isoform.pl $gff3_file";
&process_cmd($cmd);

# retrieve the best entries:
my $best_entries = "best_candidates.gff3";
$cmd = "$FindBin::Bin/../misc_utilities/gene_list_to_gff.pl $acc_file $gff3_file.inx > $best_entries";
&process_cmd($cmd);

# make a peptide file:
my $best_pep_file = $best_entries;
$best_pep_file =~ s/\.gff3$/\.pep/;
$cmd = "$FindBin::Bin/../misc_utilities/gff3_file_to_proteins.pl $best_entries $transcripts_file > $best_pep_file";
&process_cmd($cmd);

# make a CDS file:
my $best_cds_file = $best_pep_file;
$best_cds_file =~ s/\.pep$/\.cds/;
$cmd = "$FindBin::Bin/../misc_utilities/gff3_file_to_proteins.pl $best_entries $transcripts_file > $best_cds_file";
&process_cmd($cmd);

## make a BED file for viewing in IGV
my $bed_file = $best_entries;
$bed_file =~ s/\.gff3$/\.bed/;
$cmd = "$FindBin::Bin/../misc_utilities/gff3_file_to_bed.pl $best_entries > $bed_file";
&process_cmd($cmd);



exit(0);




####
sub process_cmd {
	my ($cmd) = @_;

	print "CMD: $cmd\n";
	my $ret = system($cmd);

	if ($ret) {
		die "Error, cmd: $cmd died with ret $ret";
	}
	
	return;

}





