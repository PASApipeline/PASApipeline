#!/usr/bin/env perl

use strict;
use warnings;

use PSL_parser;

my $usage = "usage: $0 blat.pslx\n\n";

my $psl_file = $ARGV[0] or die $usage;

main: {

	
	my $psl_reader = new PSL_parser($psl_file);

	while (my $entry = $psl_reader->get_next()) {

		my $mismatch_count = $entry->get_mismatch_count();
		
		my $match_count = $entry->get_match_count();

		my $gap_count = $entry->get_Q_gap_bases() + $entry->get_T_gap_bases();


		if ($mismatch_count + $gap_count) {

			my $line = $entry->get_line();

			print "$line";
		}
	}


	exit(0);

}


		
