#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use PSL_parser;

my $usage = "usage: $0 blat.psl [min_per_id=0]\n\n";

my $psl_file = $ARGV[0] or die $usage;
my $min_per_id = $ARGV[1] || 0;

main: {
	
	my @entries;
	
	my $prev_acc = "";
	
	my $psl_parser = new PSL_parser($psl_file);
	
	while (my $entry = $psl_parser->get_next()) {

		if ($min_per_id && $entry->get_per_id() < $min_per_id) {
			next;
		}
		
		my $acc = $entry->get_Q_name();
		
		if ($acc ne $prev_acc) {
			&report_top_hit(@entries) if @entries;
			@entries = ();
		}
		
		$prev_acc = $acc;
		
		my $matches = $entry->get_match_count();
		my $mismatches = $entry->get_mismatch_count();
		my $total_insertions = $entry->get_Q_gap_bases() + $entry->get_T_gap_bases();
		
		my $score = (4 * $matches) - (5 * $mismatches) - log($total_insertions+1);
		
		$entry->{__score} = $score;
		
		push (@entries, $entry);
	
	}
	
	if (@entries) {
		&report_top_hit(@entries);
	}
	
	exit(0);
}

####
sub report_top_hit {
	my (@entries) = @_;

	@entries = reverse sort {$a->{__score}<=>$b->{__score}} @entries;

	my @top_entries = shift @entries;

	while (@entries) {
		my $entry = shift @entries;
		my ($lend, $rend) = $entry->get_Q_span();
		
		my $overlaps = 0;
		foreach my $prev_entry (@top_entries) {
			my ($prev_lend, $prev_rend) = $prev_entry->get_Q_span();
			if ($prev_lend < $rend && $prev_rend > $lend) {
				$overlaps = 1;
				last;
			}
		}
		unless ($overlaps) {
			push (@top_entries, $entry);
		}
	}

	foreach my $entry (@top_entries) {

		print $entry->get_line();
	}

	
	return;
}


		
	
