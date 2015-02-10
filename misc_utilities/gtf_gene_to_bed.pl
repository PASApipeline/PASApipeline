#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Gene_obj_indexer;
use Fasta_reader;
use GFF3_utils;
use GTF_utils;
use Carp;


my $usage = "usage: $0 file.GTF > file.GFF3\n\n";

my $gtf_file = $ARGV[0] or die $usage;

main: {
            

    my $gene_obj_indexer = {};

    print STDERR "-parsing GTF file: $gtf_file\n";
    my $asmbl_id_to_gene_list_href = &GTF_utils::index_GTF_gene_objs_from_GTF($gtf_file, $gene_obj_indexer);

    foreach my $asmbl_id (sort keys %$asmbl_id_to_gene_list_href) {
        
        my @gene_ids = @{$asmbl_id_to_gene_list_href->{$asmbl_id}};
        
        #print "ASMBL: $asmbl_id, gene_ids: @gene_ids\n";
        
        foreach my $gene_id (@gene_ids) {
            
            my $gene_obj_ref = $gene_obj_indexer->{$gene_id} or die "Error, no gene obj for $gene_id";
            
            my $bed_text = $gene_obj_ref->to_BED_format();
            
            print "$bed_text";
            
        }
    }


    exit(0);
}

