#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use PSL_parser;

my $usage = "\nusage: $0 file.psl(x) MIN_PER_ID=95 MAX_PERCENT_GAP=5 STRAND=[PLUS(default),MINUS,BOTH]\n\n";

my $psl_file = $ARGV[0] or die $usage;
my $min_per_id = $ARGV[1] || 95;
my $max_percent_gap = $ARGV[2] || 5;
my $restrict_strand = $ARGV[3] || "PLUS";

main: {

	my $psl_parser = new PSL_parser($psl_file);
	
	while (my $pe = $psl_parser->get_next()) {

		my $strand = $pe->get_strand();
		if ($restrict_strand eq "PLUS" && $strand eq "-"
			||
			$restrict_strand eq "MINUS" && $strand eq '+') {
			
			next;
			
		}
		
		my $matches = $pe->get_match_count();
		my $mismatches = $pe->get_mismatch_count();

		my $gaps = $pe->get_T_gap_bases() + $pe->get_Q_gap_bases();
		
		my $per_id = $pe->get_per_id();

		my $percent_gaps = $gaps/ ($matches + $mismatches) * 100;

		#print $pe->toString() . "percent_gaps: $percent_gaps\n";
		
		#next;

		if ($per_id < $min_per_id
			||
			($percent_gaps > $max_percent_gap) ) {
			
			next;
		}

		else {

			print $pe->get_line();
		}
	}

	exit(0);
}


			
