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

my $usage = "\n\nusage: $0 gff3_file contig lend rend top_strand_Y bottom_strand_Y\n\n";

my $gff3_file = $ARGV[0] or die $usage;
my $contig = $ARGV[1] or die $usage;
my $lend = $ARGV[2] or die $usage;
my $rend = $ARGV[3] or die $usage;
my $top_strand_Y = $ARGV[4] or die $usage;
my $bottom_strand_Y = $ARGV[5] or die $usage;

my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);


my @gene_ids = @{$contig_to_gene_list_href->{$contig}};

foreach my $gene_id (@gene_ids) {
	my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};

	my ($gene_lend, $gene_rend) = sort {$a<=>$b} $gene_obj_ref->get_coords();
	
	if ($gene_lend >= $lend && $gene_rend <= $rend) {
		
		my $gene_orient = $gene_obj_ref->get_orientation();
		my $geneYpos = ($gene_orient eq '+') ? 0.75*$top_strand_Y : 0.75*$bottom_strand_Y;
		my $exonYpos = ($gene_orient eq '+') ? $top_strand_Y : $bottom_strand_Y;
		$gene_lend -= $lend;
		$gene_rend -= $lend;
		
		print "$gene_lend\t$geneYpos\n$gene_rend\t$geneYpos\n\n";
		
		my @exons = $gene_obj_ref->get_exons();
		foreach my $exon (@exons) {
			my ($exon_lend, $exon_rend) = sort {$a<=>$b} $exon->get_coords();
			
			$exon_lend -= $lend;
			$exon_rend -= $lend;

			print "$exon_lend\t$exonYpos\n$exon_rend\t$exonYpos\n\n";
		}
	}
	
	
}


exit(0);
