#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 blat.psl\n\n";

my $psl_file = $ARGV[0] or die $usage;


main: {
	
	my @entries;
	
	my $prev_acc = "";
	
	open (my $fh, $psl_file) or die "Error, cannot open file $psl_file";
	while (<$fh>) {
		unless (/^\d+\s/) {
			next;
		}
		
		my $line = $_;
		chomp;
		
		my @x = split(/\t/);
		my $acc = $x[9];
		if ($acc ne $prev_acc) {
			&report_top_hit(@entries) if @entries;
			@entries = ();
		}
		
		$prev_acc = $acc;
		
		my $matches = $x[0];
		my $mismatches = $x[1];
		my $total_insertions = $x[5]+$x[7];
		
		my $score = (4 * $matches) - (5 * $mismatches) - log($total_insertions+1);
		
		my $struct = { line => $line,
					   score => $score,
				   };
		
		push (@entries, $struct);

		#print $acc . "\n";
	}
	
	if (@entries) {
		&report_top_hit(@entries);
	}
	
	exit(0);
}

####
sub report_top_hit {
	my (@entries) = @_;

	@entries = reverse sort {$a->{score}<=>$b->{score}} @entries;

	print $entries[0]->{line};

	return;
}


		
	
