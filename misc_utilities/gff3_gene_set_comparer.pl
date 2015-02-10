#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use CdbTools;
use GFF3_utils;
use Carp;
use Nuc_translator;

my $usage = "\n\nusage: $0 before_genes.gff3 after_genes.gff3\n\n";

my $before_gff3_file = $ARGV[0] or die $usage;
my $after_gff3_file = $ARGV[1] or die $usage;


main: {
	
	my $before_gene_obj_indexer_href = {};
	my $before_contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($before_gff3_file, $before_gene_obj_indexer_href);
	my %before_isoform_genes = &convert_to_isoform_map($before_gene_obj_indexer_href);
	
	my $after_gene_obj_indexer_href = {};
	my $after_contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($after_gff3_file, $after_gene_obj_indexer_href);
	my %after_isoform_genes = &convert_to_isoform_map($after_gene_obj_indexer_href);
	
	
	foreach my $isoform_id (keys %before_isoform_genes) {

		my $before_isoform_gene = $before_isoform_genes{$isoform_id};

		my $after_isoform_gene = $after_isoform_genes{$isoform_id};
		
		if ($after_isoform_gene) {

			&compare_before_to_after($before_isoform_gene, $after_isoform_gene);
			
		}
	}
	

	exit(0);
}


####
sub compare_before_to_after {
	my ($before_gene, $after_gene) = @_;


	my $gene_id = $before_gene->{TU_feat_name};
	my $isoform_id = $before_gene->{Model_feat_name};
	
	## examine UTRs.

	my @before_5prime_utr_coords = $before_gene->get_5prime_UTR_coords();
	my @before_3prime_utr_coords = $before_gene->get_3prime_UTR_coords();


	my @after_5prime_utr_coords = $after_gene->get_5prime_UTR_coords();
	my @after_3prime_utr_coords = $after_gene->get_3prime_UTR_coords();

	my $before_5prime_utr_length = &sum_coords(@before_5prime_utr_coords);
	my $before_3prime_utr_length = &sum_coords(@before_3prime_utr_coords);

	my $after_5prime_utr_length = &sum_coords(@after_5prime_utr_coords);
	my $after_3prime_utr_length = &sum_coords(@after_3prime_utr_coords);

	if ($before_5prime_utr_length) {

		if ($after_5prime_utr_length > $before_5prime_utr_length) {
			
			my $extension = $after_5prime_utr_length - $before_5prime_utr_length;
			print join("\t", $gene_id, $isoform_id, "5primeUTR_extension", $extension) . "\n";
		}
		
	}
	elsif ($after_5prime_utr_length) {
		print join("\t", $gene_id, $isoform_id, "5primeUTR_addition", $after_5prime_utr_length) . "\n";
	}

	if ($before_3prime_utr_length) {

		if ($after_3prime_utr_length > $before_3prime_utr_length) {
			
			my $extension = $after_3prime_utr_length - $before_3prime_utr_length;

			print join("\t", $gene_id, $isoform_id, "3primeUTR_extension", $extension) . "\n";
		}
	}
	elsif ($after_3prime_utr_length) {
		print join("\t", $gene_id, $isoform_id, "3primeUTR_addition", $after_3prime_utr_length) . "\n";
	}
	
	
	return;
}


####
sub sum_coords {
	my (@coords) = @_;

	my $sum = 0;
	foreach my $coordset (@coords) {
		my ($lend, $rend) = sort {$a<=>$b} @$coordset;
		$sum += $rend - $lend + 1;
	}

	return($sum);
}


####
sub convert_to_isoform_map {
	my ($gene_map_href) = @_;

	my %isoform_map;

	foreach my $gene_id (keys %$gene_map_href) {
		
		my $gene_obj = $gene_map_href->{$gene_id};
		
		foreach my $isoform ($gene_obj, $gene_obj->get_additional_isoforms()) {

			my $isoform_id = $isoform->{Model_feat_name};

			$isoform_map{$isoform_id} = $isoform;
		}
	}


	return(%isoform_map);
}

