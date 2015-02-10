#!/usr/bin/env perl

use strict;
use warnings;


my @lengths;
while (<>) {
	chomp;
	if (/\d/) {
		push (@lengths, $_);
	}
}

@lengths = reverse sort {$a<=>$b} @lengths;

my $total_length = 0;
foreach my $len (@lengths) {
	$total_length += $len;
}

my $cummulative_lengths = 0;
foreach my $len (@lengths) {
	$cummulative_lengths += $len;
	if ($cummulative_lengths >= $total_length / 2) {
		print "N50 length: $len\n";
		last;
	}
}


exit(0);

