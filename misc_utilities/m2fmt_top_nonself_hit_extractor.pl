#!/usr/bin/env perl

use strict;
use warnings;

my $best_hit_row = "";
my $best_score = 0;

while (<STDIN>) {
	my @x = split (/\t/);
	if ($x[0] ne $x[1]) {
		if ($x[4] > $best_score) {
			$best_hit_row = $_;
			$best_score = $x[4];
		}
	}
}

if ($best_hit_row) {
	print $best_hit_row;
}


exit(0);

