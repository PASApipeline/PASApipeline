#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use GFF3_utils;
use Gene_obj;

my $usage = "usage: $0 annotation.gff3\n\n";

my $annot_gff3 = $ARGV[0] or die $usage;



main: {
	
	my $gene_obj_indexer_href = {};
	
	my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($annot_gff3, $gene_obj_indexer_href);

	foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
		
		my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};

		foreach my $gene_id (@gene_ids) {

			my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};

            print $gene_obj_ref->to_BED_format();

		}
	}


	exit(0);
	
}



