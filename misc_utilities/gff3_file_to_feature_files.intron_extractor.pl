#!/usr/bin/env perl

use strict;
use warnings;


while (<>) {
	chomp; 
	my @x = split(/\t/);
	my $contig = $x[0];
	my $lend = $x[3];
	my $rend = $x[4];
	my $orient = $x[5];
	my $seq = $x[8];

	my $upstream = substr($seq, 0, 2);
	my $downstream = substr($seq, length($seq)-2, 2);

	print join("_", $contig, $lend, $rend, $orient) . "\t$upstream....$downstream" . "\n";
}


exit(0);
