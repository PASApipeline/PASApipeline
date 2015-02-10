#!/usr/bin/env perl

use strict;
use warnings;

use Fasta_reader;

my $usage = "usage: $0 fasta_file output_prefix gapLength\n\n";

my $fasta_file = $ARGV[0] or die $usage;

my $output_prefix = $ARGV[1] or die $usage;

my $gapLength = $ARGV[2] or die $usage;

my $fasta_reader = new Fasta_reader($fasta_file);

my $bigseq = "";

open (my $ofasta_fh, ">$output_prefix.fasta") or die $!;
open (my $obed_fh, ">$output_prefix.bed") or die $!;

while (my $seq_obj = $fasta_reader->next()) {
	
	my $accession = $seq_obj->get_accession();
	my $seq = $seq_obj->get_sequence();
	
	print STDERR "\racc: $accession     ";
	
	my $before = length($bigseq) + 1;

	$bigseq .= $seq;

	my $after = length($bigseq);

	print $obed_fh join("\t", "$output_prefix", $before -1, $after -1, $accession, 0, '+',
			   $before-1, $after-1, 0, 1, $after - $before + 1, 0) . "\n";
	
	$bigseq .= 'n' x $gapLength;
	
}

print $ofasta_fh ">$output_prefix\n$bigseq\n";

print STDERR "\n\nDone.\n\n";

exit(0);

