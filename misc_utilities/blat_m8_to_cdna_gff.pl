#!/usr/bin/env perl

use strict;
use warnings;


my %key_to_segments;

while (<>) {
	chomp;
	my @x = split(/\t/);

	my $query_acc = $x[0];
	my $genome_acc = $x[1];
	
	my ($query_lend, $query_rend) = ($x[6], $x[7]);
	my ($genome_lend, $genome_rend) = ($x[8], $x[9]);

	my $Evalue = $x[10];
	
	my $bit_score = $x[11];

	
	
