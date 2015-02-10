#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;


my $usage = "\n\nusage: $0 fasta_filename seqsPerChunk core_filename\n\n"; 

my $fasta_filename = $ARGV[0] or die $usage;
my $seqsPerChunk = $ARGV[1] or die $usage;
my $core_filename = $ARGV[2] or die $usage;

main:{
	

	my $seq_counter = 0;
	my $ofh;
	
	my $fasta_reader = new Fasta_reader($fasta_filename);
	while (my $seq_obj = $fasta_reader->next()) {
		
		my $fasta_format = $seq_obj->get_FASTA_format();
		$seq_counter++;
		if ($seq_counter % $seqsPerChunk == 1) {
			close $ofh if $ofh;
			open ($ofh, ">$core_filename.$seq_counter.fasta") or die $!;
			print STDERR "-writing file: $core_filename.$seq_counter.fasta\n";
		}
		print $ofh $fasta_format . "\n";
	}
	close $ofh if $ofh;


	print "Done.\n\n";
}
		

	   
