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


my $usage = "usage: $0 file.GFF3 > file.single_iso.GFF3\n\n";

my $gff3_file = $ARGV[0] or die $usage;


main: {
    my $gene_obj_indexer = {};
    
    my $asmbl_id_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer);

    foreach my $asmbl_id (sort keys %$asmbl_id_to_gene_list_href) {
        
        my @gene_ids = @{$asmbl_id_to_gene_list_href->{$asmbl_id}};
        
        #print "ASMBL: $asmbl_id, gene_ids: @gene_ids\n";
        
        foreach my $gene_id (@gene_ids) {
            
            ## note models of isoforms are bundled into the same gene object.
            my $gene_obj_ref = $gene_obj_indexer->{$gene_id} or die "Error, no gene retrieved based on id: $gene_id";
            
            my @gene_and_cds_lengths;
            
            foreach my $gene_obj ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
                
                $gene_obj->delete_isoforms(); # unbundle the model object!

                my $cds_len = $gene_obj->get_CDS_length();

                push (@gene_and_cds_lengths, [$cds_len, $gene_obj]);
                
            }

            @gene_and_cds_lengths = reverse sort {$a->[0] <=> $b->[0]} @gene_and_cds_lengths;

            print $gene_and_cds_lengths[0]->[1]->to_GFF3_format() . "\n";
            
        }
    }

    
    exit(0);
}

