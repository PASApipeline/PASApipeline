#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Nuc_translator;
use Fasta_reader;

my $usage = "usage: $0 fasta_file lend rend\n";

my $fasta_file = $ARGV[0] or die $usage;
my $lend = $ARGV[1] or die $usage;
my $rend = $ARGV[2] or die $usage;

unless ($lend < $rend) { die $usage;}


main: {
	
	my $fasta_reader = new Fasta_reader($fasta_file);
	
	my $seq_obj = $fasta_reader->next();

	my $header = $seq_obj->get_header();
	my $sequence = lc $seq_obj->get_sequence();
	
	my @seq_chars = split (//, $sequence);
	for (my $i = $lend; $i <= $rend; $i++) { 
		
		$seq_chars[$i-1] = uc $seq_chars[$i-1];
	}


	$sequence = join ("", @seq_chars);

	$sequence =~ s/(\S{60})/$1\n/g;

	print ">$header\n$sequence\n";

	exit(0);
}
	
