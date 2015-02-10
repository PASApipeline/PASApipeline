#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "$0 m2fmt_file per_ID covA covB\n\n";

my $m2fmt_file = $ARGV[0] or die $usage;
my $per_ID = $ARGV[1] or die $usage;
my $covA = $ARGV[2] or die $usage;
my $covB = $ARGV[3] or die $usage;

main: {
	open (my $fh, $m2fmt_file) or die "Error, cannot open file $m2fmt_file";
	while (<$fh>) {
		chomp;
		my @x = split (/\t/);
		if ($x[10] >= $per_ID && $x[23] >= $covA && $x[25] >= $covB) {
			print "$_\n";
		}
	}
	close $fh;
	
	exit(0);
}



