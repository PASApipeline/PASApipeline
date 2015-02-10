#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 file.psl [min_query_length=0]\n\n";

my $blat_file = $ARGV[0] or die $usage;
my $min_query_length = $ARGV[1] || 0;


my $min_per_id = 98;


my $total_match = 0;
my $total_mismatch = 0;
my $total_insertions = 0;
my $total_deletions = 0;

open (my $fh, $blat_file) or die "Error, cannot open file $blat_file";
while (<$fh>) {
	unless (/\w/) { next; }
	
	unless (/^\d+\s/) { next; }
	
	
	chomp;
	my @x = split(/\t/);
	
	my $query_len = $x[10];
	unless ($query_len >= $min_query_length) { 
		next;
	}
	
	my $match_count = $x[0];
	my $mismatch_count = $x[1];
	my $asmbl_gaps = $x[5];
	my $gene_gaps = $x[7];
	
	if ($match_count / ($match_count + $mismatch_count) * 100 < $min_per_id) { 
		next;
	}
	

	$total_match += $match_count;
	$total_mismatch += $mismatch_count;
	$total_deletions += $gene_gaps;
	$total_insertions += $asmbl_gaps;
	
	
}

print "Match: $total_match\n"
	. "Mismatch: $total_mismatch\n"
	. "Insertions: $total_insertions\n"
	. "Deletions: $total_deletions\n";


print "\n% mismatches: " . sprintf("%.6f", $total_mismatch/($total_mismatch+$total_match)*100) 
	. "\t%insertions: " . sprintf("%.6f", $total_insertions/($total_mismatch+$total_match+$total_insertions)*100)
	. "\n\n";


exit(0);

