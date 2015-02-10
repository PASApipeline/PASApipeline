#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Overlap_piler;

my $usage = "\nusage: $0 contig_coords.file\n\nFormat should be\nContig(tab)end5(tab)end3\n...\n\n";

my $file = $ARGV[0] or die $usage;

main: {
	
	my %contig_to_coordpairs;
	
	open (my $fh, $file) or die "Error, cannot open file $file";
	while (<$fh>) {
		chomp;
		my ($contig, $lend, $rend) = split (/\t/);
		unless ($lend =~ /^\d+$/ && $rend =~ /^\d+$/) {
			die "Error, couldn't parse line $_";
		}
		
		push (@{$contig_to_coordpairs{$contig}}, [$lend, $rend]);
	}
	close $fh;

	foreach my $contig (sort keys %contig_to_coordpairs) {
		my $coords_list_aref = $contig_to_coordpairs{$contig};

		my @clusters = &Overlap_piler::simple_coordsets_collapser(@$coords_list_aref);
		foreach my $cluster (@clusters) {
			my ($lend, $rend) = @$cluster;
			print "$contig\t$lend\t$rend\n";
		}
	}

	exit(0);
}



