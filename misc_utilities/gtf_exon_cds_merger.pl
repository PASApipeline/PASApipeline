#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Gene_obj_indexer;
use CdbTools;
use GFF3_utils;
use GTF_utils;
use Carp;
use Overlap_piler;
use Data::Dumper;

my $usage = "usage: $0 file.GTF genome.fasta > file.GTF\n\n";

my $gtf_file = $ARGV[0] or die $usage;
my $genome_fasta = $ARGV[1] or die $usage;


main: {
	

	my $asmbl_id_to_genes_and_components_href = &parse_exons_and_cds($gtf_file);

	foreach my $asmbl_id (keys %$asmbl_id_to_genes_and_components_href) {
		
		my $genome_seq = cdbyank_linear($asmbl_id, $genome_fasta);
		
		my $genes_href = $asmbl_id_to_genes_and_components_href->{$asmbl_id};

		foreach my $gene_id (keys %$genes_href) {

			my $orient = $genes_href->{$gene_id}->{orient};

			my $transcripts_href = $genes_href->{$gene_id}->{structure};

			## get all the end5, end3 values.
			## cluster them by overlap, and recreate a gene object
			
			foreach my $transcript_id (keys %$transcripts_href) {
				
				my $coords_aref = $transcripts_href->{$transcript_id};

				my @merged_coords = &Overlap_piler::simple_coordsets_collapser(@$coords_aref);
				
			
				## convert to end5, end3 values:
				my %coords;
				foreach my $coordset (@merged_coords) {
					my ($lend, $rend) = @$coordset;
				    my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);
					
					$coords{$end5} = $end3;
				}


				## convert to gene object and output GTF
				my $gene_obj = new Gene_obj();
				$gene_obj->populate_gene_object(\%coords, \%coords);
				$gene_obj->{TU_feat_name} = $gene_id;
				$gene_obj->{Model_feat_name} = $transcript_id;
				$gene_obj->{asmbl_id} = $asmbl_id;
				$gene_obj->{com_name} = "";
				print $gene_obj->to_GFF3_format() . "\n";
				
				
			}
		}
	}
	
	

	exit(0);
}
			



####
sub parse_exons_and_cds {
	my ($gtf_file) = @_;

	my %data;

	open (my $fh, $gtf_file) or die "error, cannot open file $gtf_file";
	while (<$fh>) {
		unless (/\w/) { next; }
		if (/^\#/) { next; }
		chomp;
		my ($contig, $type, $feat, $lend, $rend, $score, $orient, $phase, $gene_info) = split (/\t/);

		
		if ($feat eq 'CDS' || $feat eq 'exon') {
			
			$gene_info =~ /gene_id \"([^\"]+)/;
			my $gene_id = $1 or die "Error, cannot parse gene ID from $gene_info";

			$gene_info =~ /transcript_id \"([^\"]+)/;
			my $transcript_id = $1 or die "Error, cannot parse transcript ID from $gene_info";

			push (@{$data{$contig}->{$gene_id}->{structure}->{$transcript_id}}, [$lend, $rend]);
			$data{$contig}->{$gene_id}->{orient} = $orient;
			
		}
		
	}
	
	return (\%data);
}

			
