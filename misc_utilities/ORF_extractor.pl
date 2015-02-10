#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling);

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;
use Longest_orf;

my $usage = <<_EOUSAGE_;

####################################################################################
#
# Required:
#
#  --fasta <filename>    fasta file
#
# Optional:
#
#  --allow_partials            allows both 5' and 3' partials (assumes --allow_3prime_partials and --allow_5prime_partials)
#  --allow_3prime_partials    
#  --allow_5prime_partials 
#
#  --forward_strand_only
#  --reverse_strand_only
#
#  --allow_non_met_starts       finds longest orf from stop+1 to stop
#
#  --min_ORF_length            default: 10 amino acids
#  --max_number_ORFs           unset, by default shows all.
#
#################################################################################

_EOUSAGE_

	;

my $help_flag;

# param opts
my ($fasta_file);
my $min_ORF_length = 10;

# flag opts
my ($allow_partials_flag, $allow_3prime_partials_flag, $allow_5prime_partials_flag,
	$forward_strand_only_flag, $reverse_strand_only_flag,
	$allow_non_met_starts_flag, $max_number_ORFs,
	);



&GetOptions ( 'h' => \$help_flag,

			  'fasta_file=s' => \$fasta_file,
			  'min_ORF_length=i' => \$min_ORF_length,


			  'allow_partials' => \$allow_partials_flag,
			  'allow_3prime_partials' => \$allow_3prime_partials_flag,
			  'allow_5prime_partials' => \$allow_5prime_partials_flag,
			  'forward_strand_only' => \$forward_strand_only_flag,
			  'reverse_strand_only' => \$reverse_strand_only_flag,
			  'allow_non_met_starts' => \$allow_non_met_starts_flag,
			  'max_number_ORFs=i' => \$max_number_ORFs,
			  
	);

unless ($fasta_file && -s $fasta_file) {
	die $usage;
}


main: {
	my $longest_orf_finder = new Longest_orf();

	if ($allow_partials_flag) {
		$longest_orf_finder->allow_partials();
	}
	else {
		if ($allow_5prime_partials_flag) {
			$longest_orf_finder->allow_5prime_partials();
		}
		if ($allow_3prime_partials_flag) {
			$longest_orf_finder->allow_3prime_partials();
		}
	}

	if ($forward_strand_only_flag) {
		$longest_orf_finder->forward_strand_only();
	}
	elsif ($reverse_strand_only_flag) {
		$longest_orf_finder->reverse_strand_only();
	}
	
	if ($allow_non_met_starts_flag) {
		$longest_orf_finder->allow_non_met_starts();
	}

	## do the work
	
	my $fasta_reader = new Fasta_reader($fasta_file);
	
	while (my $seq_obj = $fasta_reader->next()) {
		
		my $accession = $seq_obj->get_accession();
		my $sequence = $seq_obj->get_sequence();
		
		my @orf_structs = $longest_orf_finder->capture_all_ORFs($sequence);

		@orf_structs = reverse sort {$a->{length}<=>$b->{length}} @orf_structs;
		
		if ($max_number_ORFs && scalar(@orf_structs) > $max_number_ORFs) {
			@orf_structs = @orf_structs[0..$max_number_ORFs-1];
		}
		
		my @results;
		
		foreach my $orf (@orf_structs) {
			my $start = $orf->{start};
			my $stop = $orf->{stop};
			my $length = $orf->{length}/3;
			my $orient = $orf->{orient};
			my $protein = $orf->{protein};
			
			if ($length >= $min_ORF_length) {
				
				push (@results, "$length:$start-$stop($orient),$protein");
			}
		}

		print "$accession\t" . join("\t", @results) . "\n";
		
		
	}
	
	

	exit(0);
}

