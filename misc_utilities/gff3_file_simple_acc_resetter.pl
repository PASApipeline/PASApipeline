#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use CdbTools;
use GFF3_utils;
use Carp;

my $usage = "\n\nusage: $0 gff3_file\n\n";

my $gff3_file = $ARGV[0] or die $usage;

my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

my $gene_counter = 0;

foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
    
    my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
    foreach my $gene_id (@gene_ids) {
        my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};

		$gene_counter++;

		my $isoform_counter = 0;
		
		foreach my $isoform ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
			
			$isoform_counter++;
			
			$isoform->{TU_feat_name} = "Gene.$gene_counter";
			$isoform->{Model_feat_name} = "mRNA.$gene_counter.$isoform_counter";
			
		}
		
		print $gene_obj_ref->to_GFF3_format();
		print "\n";
		
	}
}



exit(0);

