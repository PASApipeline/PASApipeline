#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 btab_files ortholog_summary\n\n";

my $btab_file = $ARGV[0] or die $usage;
my $ortholog_summary = $ARGV[1] or die $usage;

my %ortholog_pairs;

{ # parse the ortholog pairs
	
	open (my $fh, $ortholog_summary) or die "error, cannot open file $ortholog_summary";
	while (<$fh>) {
		chomp;
		my @x = split (/\t/);
		my ($orthoA, $orthoB) = ($x[1], $x[2]);
		
		my $ortho_token = join ("$;", $orthoA, $orthoB);

		$ortholog_pairs{$ortho_token} = 1;
	}
	close $fh;
}

{ # parse the btab
	
	open (my $fh, $btab_file) or die "error, cannot open file $btab_file";
	while (<$fh>) {
		my $line = $_;
		
		chomp;
		my @x = split (/\t/);
		my ($accA, $accB) = ($x[0], $x[5]);

		my $token1 = join ("$;", $accA, $accB);
		my $token2 = join ("$;", $accB, $accA);

		if ($ortholog_pairs{$token1} || $ortholog_pairs{$token2}) {
			print $line;
		}
	}
	close $fh;
}


exit(0);



	
	
	
		
