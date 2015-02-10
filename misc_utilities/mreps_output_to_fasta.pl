#!/usr/bin/env perl

use strict;
use warnings;

while (<>) {
	chomp;
	s/^\s+//;
	my ($contig, $lend, $dash, $rend, $colon, $length, $period, $copies, $variance, @repeats) = split (/\s+/);

  
	my %repeat_units;
	foreach  my $repeat (@repeats) {
		$repeat_units{$repeat}++;
	}

	my @reps = sort {$repeat_units{$a}<=>$repeat_units{$b}} keys %repeat_units;
	
	my $representative_repeat_unit = pop @reps;

	my $repeat_string = join ("", @repeats);
	
	print ">$contig-$lend-$rend $period $representative_repeat_unit $variance\n$repeat_string\n";

}


exit(0);
