#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;

my $usage = "usage: $0 translationsFasta min_length\n\n";

my $fasta_file = $ARGV[0] or die $usage;
my $min_length = $ARGV[1] or die $usage;

main: {

	my $fasta_reader = new Fasta_reader($fasta_file);
	
	while (my $seq_obj = $fasta_reader->next()) {
		
		my $acc = $seq_obj->get_accession();
		my $sequence = uc $seq_obj->get_sequence();

		my @chars = split (//, $sequence);
		
		for (my $i = 0; $i <= $#chars - $min_length; $i++) {
			
			if ($i==0 || $chars[$i] eq 'M') {
				# walk to a stop codon or end of the sequence
				my $j = $i;
				
				my @chars_in_range;
				while ($j <= $#chars&& $chars[$j] ne '*') { 
					push (@chars_in_range, $chars[$j]);
					$j++;
				}
				if (scalar @chars_in_range >= $min_length) {
					
					my $orf = join ("", @chars_in_range);
					print ">$acc-$i\n$orf\n";
				}
			}
		}
	}
	exit(0);

}



