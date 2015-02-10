#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;
use CdbTools;

my $usage = "\nusage: $0 nuc.fasta max_percent_repeat\n\n";

my $nuc_fasta_file = $ARGV[0] or die $usage;
my $max_percent_repeat = $ARGV[1] or die $usage;



main: {


	my @accessions_remaining;

	my $tmpfile = "$$.nseg";

	my $cmd = "nseg $nuc_fasta_file -x > $tmpfile";
	my $ret = system $cmd;
	if ($ret) {
		die "Error, cmd: $cmd died with ret $ret";
	}
	
	my $fasta_reader = new Fasta_reader($tmpfile);

	while (my $seq_obj = $fasta_reader->next()) {
		
		my $accession = $seq_obj->get_accession();
		
		my $sequence = $seq_obj->get_sequence();
		my $seq_length = length($sequence);
		
		my $num_N = 0;
		while ($sequence =~ /n/ig) { $num_N++;}
		
		my $percent_N = $num_N / $seq_length * 100;

		if ($percent_N <= $max_percent_repeat) {
			push (@accessions_remaining, $accession);
		}

	}
	

	foreach my $acc (@accessions_remaining) {
		
		my $fasta_entry = cdbyank($acc, $nuc_fasta_file);
		
		print $fasta_entry;
	}
	
	exit(0);
}

