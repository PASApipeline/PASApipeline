#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 neighbors.out.txt\n\n";

my $neighbors_file = $ARGV[0] or die $usage;


my %best_hits;

open (my $fh, $neighbors_file) or die "Error, cannot open file $neighbors_file";
while (<$fh>) {
	chomp;
	
	my ($accA, $accB, $hsp_len_A, $hsp_len_B, $blast_score, $per_id) = split(/\t/);
	
	if ($accA eq $accB) { next; }
	
	foreach my $acc ($accA, $accB) {
		
		my $struct = $best_hits{$acc};
		
		if ( (! $struct) || $struct->{score} < $blast_score) {
		
			$best_hits{$acc} = { accA => $accA,
								 accB => $accB,
								 score => $blast_score,
								 per_id => $per_id,
							 };

		}

	}
}

close $fh;

foreach my $struct (values %best_hits) {
	
	print join("\t", 
			   $struct->{accA},
			   $struct->{accB},
			   $struct->{score},
			   $struct->{per_id}) . "\n";
}


exit(0);


		
