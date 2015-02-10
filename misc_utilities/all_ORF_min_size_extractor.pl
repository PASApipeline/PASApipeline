#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;

my $usage = "usage: $0 translationsFasta min_length\n\n";

my $fasta_file = $ARGV[0] or die $usage;
my $min_length = $ARGV[1] or die $usage;

$min_length--;  #adding a non-stop in the regex

main: {
	my $fasta_reader = new Fasta_reader($fasta_file);
	
	while (my $seq_obj = $fasta_reader->next()) {
		
		my $acc = $seq_obj->get_accession();
		my $sequence = $seq_obj->get_sequence();

		while ($sequence =~ /(M[^\*X]{$min_length}[^\*X]*)/g) {
			my $start = $-[0];
			my $subseq = $1;
			
			print ">$acc-$start\n$subseq\n";
		}
	}

	exit(0);

}



