#!/usr/bin/env perl

use strict;
use warnings;

use Fasta_reader;
use Nuc_translator;


my $usage = "usage: $0 nucFasta\n\n";

my $fasta_file = $ARGV[0] or die $usage;


main: {
	my $fasta_reader = new Fasta_reader($fasta_file);
	while (my $seq_obj = $fasta_reader->next()) {
		
		my $accession = $seq_obj->get_accession();
		my $sequence = $seq_obj->get_sequence();

		foreach my $frame (1,2,3) {
			my $translation = &translate_sequence($sequence, $frame);
			
			$translation =~ s/(\S{60})/$1\n/g;
			chomp $translation;
			
			print ">$accession^Fp$frame\n$translation\n";
		}


		$sequence = &reverse_complement($sequence);
		
		foreach my $frame (1,2,3) {
			my $translation = &translate_sequence($sequence, $frame);
			
			$translation =~ s/(\S{60})/$1\n/g;
			chomp $translation;
			
			print ">$accession^Fm$frame\n$translation\n";
		}
		
	}

	exit(0);
}



