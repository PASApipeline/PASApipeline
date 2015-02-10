#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use MALIGN::Alignment_parser;
use MALIGN::Alignment;
use MALIGN::Aligned_sequence;

my $usage = "usage: $0 alignment_file format\n\n"
	. " allowable formats include: stockholm, phylip\n\n";


my $alignment_file = $ARGV[0] or die $usage;
my $format = $ARGV[1] or die $usage;

main: {
	my $alignment_obj = MALIGN::Alignment_parser::parse_alignment($alignment_file, $format);
	
	my @aligned_seqs = $alignment_obj->get_aligned_seqs();

	## compare percent identities
	
	for (my $i = 0; $i < $#aligned_seqs; $i++) {

		my $aligned_i = $aligned_seqs[$i];
		
		for (my $j = $i + 1; $j <= $#aligned_seqs; $j++) {
			
			my $aligned_j = $aligned_seqs[$j];

			my $percent_identity = $aligned_i->percent_identity_to($aligned_j);

			print $aligned_i->get_accession() . "\t" . $aligned_j->get_accession() . "\t" . sprintf("%.2f", $percent_identity) . "\n";
			
		}
	}


	exit(0);
}



		
		
