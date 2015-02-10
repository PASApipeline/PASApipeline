#!/usr/bin/env perl

use strict;
use warnings;


my $usage = "usage: $0 gff3_file\n\n";

my $gff3 = $ARGV[0] or die $usage;

main: {
	
	my %contig_to_gene_list;
	
	my $gene_counter = 0;

	open (my $fh, $gff3) or die "Error, cannot open file $gff3";
	while (<$fh>) {
		chomp;
		unless (/\w/) { next; }
		my @x = split (/\t/);
		if ($x[2] eq 'gene') {
	
			$gene_counter++;
			
			my ($scaff, $lend, $rend, $name) = ($x[0], $x[3], $x[4], $x[8]);

			push (@{$contig_to_gene_list{$scaff}}, { name => $name,
													 lend => $lend,
													 rend => $rend,
													 midpt => ($lend + $rend) / 2,
												 } );
		}
	}
	close $fh;

	
	my @intergenics;
	
	foreach my $scaffold (keys %contig_to_gene_list) {
		
		my @genes = sort {$a->{midpt}<=>$b->{midpt}} @{$contig_to_gene_list{$scaffold}};

		my $prev_gene = shift @genes;
		while (@genes) {
			my $next_gene = shift @genes;

			my $prev_acc = $prev_gene->{name};
			my $next_acc = $next_gene->{name};

			my $delta = $next_gene->{midpt} - $prev_gene->{midpt};
			# print "$scaffold\t$prev_acc\t$next_acc\t$delta\n";
			if ($delta > 0) {
				push (@intergenics, { dist => $delta,
									  accA => $prev_acc,
									  accB => $next_acc,
									  
								  } );
			}

			$prev_gene = $next_gene;
		}

	}

	@intergenics = sort {$a->{dist}<=>$b->{dist}} @intergenics;

	my $intergene_counter = 0;
	my %seen;
	my $sum_interdist = 0;
	
	my $prev_percent = 0;

	foreach my $intergenic (@intergenics) {
		my ($dist, $accA, $accB) = ($intergenic->{dist}, $intergenic->{accA}, $intergenic->{accB});
		
		$sum_interdist += $dist;
		unless ($seen{$accA}) {
			$seen{$accA} = 1;
			$intergene_counter++;
		}
		unless ($seen{$accB}) {
			$seen{$accB} = 1;
			$intergene_counter++;
		}

		
		my $percent = sprintf ("%.2f", $intergene_counter / $gene_counter * 100);
		if (int($percent) > $prev_percent && int($percent) % 5 == 0) {
			print "$percent\t$intergene_counter\t$sum_interdist\n";
			$prev_percent = int($percent);
		}
	}
	
	exit(0);
}
	
	
	
		
