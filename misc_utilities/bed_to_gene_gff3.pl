#!/usr/bin/env perl

use strict;
use warnings;


use lib ($ENV{EUK_MODULES});
use Gene_obj;
use Fasta_reader;

use Data::Dumper;

my $usage = "usage: $0 annots.bed\n\n";

my $annots_bed = $ARGV[0] or die $usage;

my $DEBUG = 0;

main: {
	
	my $counter = 0;

	open (my $fh, $annots_bed) or die "Error, cannot open file $annots_bed";
	while (<$fh>) {
		print if $DEBUG;
		
		my @x = split(/\t/);

		my $scaff = $x[0];
		my $gene_lend = $x[1] + 1;
		my $gene_rend = $x[2];

		my $com_name = $x[3];

		my $score = $x[4];
		my $orient = $x[5];
		
		if ($orient eq '*') {
			$orient = '+';
		}
		

		my $coding_lend = $x[6] + 1;
		my $coding_rend = $x[7];

		my $rgb_color = $x[8];

		my $num_exons = $x[9];

		my $lengths_text = $x[10];
		my $exon_relative_starts_text = $x[11];

		my @lengths = split(/,/, $lengths_text);
		my @exon_relative_starts = split(/,/, $exon_relative_starts_text);

		my @exons;

		my $sum_len = 0;

		while (@lengths) {
			my $len = shift @lengths;
			my $start = shift @exon_relative_starts;
			
			my $exon_lend = $gene_lend + $start;
			my $exon_rend = $exon_lend + $len - 1;
			

			print "Len: $len, start=$start   ====>  $exon_lend - $exon_rend\n" if $DEBUG;

			push (@exons, [$exon_lend, $exon_rend]);
			

			$sum_len += $len;
		}

		if ($sum_len < 3) {
			print STDERR "Ignoring entry: $_, since ultra short\n";
			next;
		}
		
		
		print "Coding: $coding_lend-$coding_rend, Exons: " . Dumper (\@exons) if $DEBUG;
		
		eval {
			
			my $gene_obj = new Gene_obj();

            if ($coding_lend == $coding_rend +1) { ## not coding
                $coding_lend = 0;
                $coding_rend = 0;
            }
            
			$gene_obj->build_gene_obj_exons_n_cds_range(\@exons, $coding_lend, $coding_rend, $orient);
			
			$gene_obj->{com_name} = $com_name;
			$gene_obj->{asmbl_id} = $scaff;
			
			$counter++;
			$gene_obj->{TU_feat_name} = "gene.$counter";
			$gene_obj->{Model_feat_name} = "model.$counter";
			
			print $gene_obj->to_GFF3_format() . "\n";
		
		};

		if ($@) {
			print "Error:\n$@\n";
		}
	}


	exit(0);
}
	
