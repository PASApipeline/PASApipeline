#!/usr/bin/env perl

use strict;
use warnings;

use Overlap_piler;

my $usage = "usage: $0 feature_file [termini_extend=0]\n\n"
	. "where feature file has format:\n\n"
	. "contig\tend5\tend3\tacc\n\n\n";

my $file = $ARGV[0] or die $usage;
my $term_extend = $ARGV[1] || 0;

main: {

	my %contig_to_features;

	open (my $fh, $file) or die "Error, cannot open file $file";
	while (<$fh>) {
		chomp;
		my ($contig, $lend, $rend, $acc) = split(/\t/);
		
		($lend, $rend) = sort {$a<=>$b} ($lend, $rend);

		push (@{$contig_to_features{$contig}}, [$acc, $lend, $rend]);

	}
	close $fh;


	foreach my $contig (keys %contig_to_features) {

		my @feats = @{$contig_to_features{$contig}};

		my $coordset_collapser = new Overlap_piler();
		
		my %acc_to_coords;

		foreach my $feat (@feats) {
			my ($acc, $lend, $rend) = @$feat;
			$coordset_collapser->add_coordSet($acc, $lend-$term_extend, $rend+$term_extend);
		
			$acc_to_coords{$acc} = [$lend, $rend];
		}

		my @clusters = $coordset_collapser->build_clusters();
		
		foreach my $cluster (@clusters) {
			
			my @accs = @$cluster;

			my @coords;
			foreach my $acc (@accs) {
				my $coords_aref = $acc_to_coords{$acc};
				push (@coords, @$coords_aref);
			}

			@coords = sort {$a<=>$b} @coords;
			
			my $span_lend = shift @coords;
			my $span_rend = pop @coords;

			print "$contig\t$span_lend\t$span_rend\t" . scalar(@accs) . "\t" . join(" ", @accs) . "\n";
			
		}

	}


	exit(0);
}

