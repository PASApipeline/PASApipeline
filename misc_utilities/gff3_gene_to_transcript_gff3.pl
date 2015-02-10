#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
#use lib ("$FindBin::Bin/../PerlLib");

use lib ($ENV{EUK_MODULES});

use Gene_obj;
use GFF3_utils;
use Carp;


my $usage = "\n\nusage: $0 gff3_file\n\n";

my $gff3_file = $ARGV[0] or die $usage;

my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);


my $trans_counter = 0;

foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
    my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
    foreach my $gene_id (@gene_ids) {
        		
		my $gene = $gene_obj_indexer_href->{$gene_id};
		
        my $id = "t.$$." . ++$trans_counter;

		my $align_gff3_text = $gene->to_alignment_GFF3_format($id, $gene_id);
			
		print $align_gff3_text . "\n";
		
	}
}


exit(0);

