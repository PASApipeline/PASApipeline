#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;

my $usage = "\nusage: $0 fasta_file\n\n";

my $fasta_file = $ARGV[0] or die $usage;
srand();

my @chars = qw (g a t c);


main: {

	

	my $fasta_reader = new Fasta_reader($fasta_file);

	while (my $seq_obj = $fasta_reader->next()) {

		my $accession = $seq_obj->get_accession();
		my $fasta_entry = $seq_obj->get_FASTA_format();
		my $sequence = $seq_obj->get_sequence();
		

		my @N_regions;
		while ($sequence =~ /(n+)/ig) {
			push (@N_regions, [$-[0], $+[0]]);
		}
		
		my @seq_chars = split (//, $sequence);
		foreach my $N_region (@N_regions) {
			my ($lend, $rend) = @$N_region;
			
			for (my $i = $lend; $i <= $rend; $i++) {
				
				my $rand_char = $chars[ int(rand(4)) ];
				$seq_chars[$i] = $rand_char;
			}
		}
		$sequence = join ("", @seq_chars);
				
		$sequence =~ s/(\S{60})/$1\n/g;
		

		open (my $fh, ">$accession.$$.seq") or die "Error, cannot write file";
		print $fh ">$accession\n$sequence\n";
		close $fh;
		
		## run mreps:
		my $cmd = "mreps -res 3 -exp 3.0 -minperiod 3 -minsize 50 -fasta $accession.$$.seq";

		open (my $cmdfh, "$cmd | ");
		while (<$cmdfh>) {
			#print;
			if (/^\s*\d+\s+\->\s+\d+/) {
				s/[<>]//g;
				print $accession . "\t$_";
			}
		}
		close $cmd;
		

		unlink("$accession.$$.seq");
		
	}

	exit(0);
}


			
