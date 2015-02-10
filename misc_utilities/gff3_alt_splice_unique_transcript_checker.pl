#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use GFF3_utils;
use Carp;
use Gene_obj_comparator;


my $usage = "\n\nusage: $0 gff3_file \n\n";

my $gff3_file = $ARGV[0] or die $usage;

my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
    
    my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
    foreach my $gene_id (@gene_ids) {
        my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
		
		my @isoforms = ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms());
		
		my $num_isoforms = scalar(@isoforms);
		if ($num_isoforms > 1) {
			
			my $gene_id = $gene_obj_ref->{TU_feat_name};

			for (my $i = 0; $i < $#isoforms; $i++) {
				for (my $j = $i + 1; $j <= $#isoforms; $j++) {
					
					my $model_i = $isoforms[$i]->{Model_feat_name};
					my $model_j = $isoforms[$j]->{Model_feat_name};
					
					&compare_genes($isoforms[$i], $isoforms[$j]);
					
					
					
					my $mRNAsame = &are_mRNA_same();
					my $CDSsame = &are_CDS_same();
					
					print "$gene_id\t$model_i\t$model_j\tmRNAsame:$mRNAsame\tCDSsame:$CDSsame";
					if ($mRNAsame && $CDSsame) {
						print "\tIDENTICAL\n";
					}
					else {
						print "\n";
					}
					
					
				}
			}
		}
	}


}

exit(0);


