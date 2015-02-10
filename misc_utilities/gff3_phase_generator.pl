#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/PerlLib", "$FindBin::Bin/../PerlLib");
use Gene_obj;
use Gene_obj_indexer;
use GFF3_utils;
use CdbTools;

my $usage = "\n\nusage: $0 gff3_file genome_db\n\n";

my $gff3_file = $ARGV[0] or die $usage;
my $fasta_db = $ARGV[1] or die $usage;


my $index_file = "$gff3_file.inx";
my $gene_obj_indexer = new Gene_obj_indexer( { "create" => $index_file } );

## associate gene identifiers with contig id's.
my $asmbl_id_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer);

foreach my $asmbl_id (keys %$asmbl_id_to_gene_list_href) {
    
    my $genome_seq = cdbyank_linear($asmbl_id, $fasta_db);
    
    my @gene_ids = @{$asmbl_id_to_gene_list_href->{$asmbl_id}};
    
    foreach my $gene_id (@gene_ids) {
        my $gene_obj_ref = $gene_obj_indexer->get_gene($gene_id);
        
        foreach my $isoform ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
            $isoform->delete_isoforms();
            
            my $isoform_id = $isoform->{Model_feat_name};
            
            
            eval {
                $isoform->set_CDS_phases (\$genome_seq);
                
                
            };
            
            if ($@) {
                print STDERR "Error encountered: $@\nCDS phases could not be established for $isoform_id\n";
            }
            
            print $isoform->to_GFF3_format() . "\n";
        }
    }
}


unlink ($index_file);

exit(0);

