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


my $usage = "usage: $0 file.GTF [genome.fasta] > file.GFF3\n\n";

my $gtf_file = $ARGV[0] or die $usage;
my $genome_fasta = $ARGV[1];


main: {
            

    my $gene_obj_indexer = {};

    print STDERR "-parsing GTF file: $gtf_file\n";
    my $asmbl_id_to_gene_list_href = &GTF_utils::index_GTF_gene_objs_from_GTF($gtf_file, $gene_obj_indexer);

    my %genome;
    if ($genome_fasta) {
        print STDERR "-parsing fasta file: $genome_fasta\n";
        my $fasta_reader = new Fasta_reader($genome_fasta);
        %genome = $fasta_reader->retrieve_all_seqs_hash();
    }
    

    foreach my $asmbl_id (sort keys %$asmbl_id_to_gene_list_href) {
        
        ## get the genome sequence
        my $genome_seq;
        if ($genome_fasta) {
            $genome_seq = $genome{$asmbl_id} or die "Error, no sequence for asmbl_id: $asmbl_id";
        }
        
        my @gene_ids = @{$asmbl_id_to_gene_list_href->{$asmbl_id}};
        
        #print "ASMBL: $asmbl_id, gene_ids: @gene_ids\n";
        
        foreach my $gene_id (@gene_ids) {
            
            my $gene_obj_ref = $gene_obj_indexer->{$gene_id} or die "Error, no gene obj for $gene_id";
            
            foreach my $gene_obj ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
                
                if (my $com_name = $gene_obj->{com_name}) {
                    $gene_obj->{TU_feat_name} = "$com_name|" . $gene_obj->{TU_feat_name};
                    $gene_obj->{Model_feat_name} = "$com_name|" . $gene_obj->{Model_feat_name};
                }
                
                                
                if ($genome_seq) {
                    eval {
                        $gene_obj->set_CDS_phases(\$genome_seq);
                    };
                    if ($@) {
                        # do it in pseudogene mode
                        print STDERR "Warning, gene obj " . $gene_obj->toString() . " may be corrupt. $@\n";
						
                    }
                }
                
                
				
				if ($genome_fasta) {
					my $cds_sequence = $gene_obj->get_CDS_sequence();
					my $protein_sequence = $gene_obj->get_protein_sequence();
					
					my $gene_feat = $gene_obj->{TU_feat_name};
					my $model_feat = $gene_obj->{Model_feat_name};
					if ($protein_sequence && $cds_sequence) {
                        print "# $model_feat $gene_feat\tPROTEIN\t$protein_sequence\n\n";
                        print "# $model_feat $gene_feat\tCDS\t$cds_sequence\n\n";
                    }
                }
            }
            
            my $gff3_text = $gene_obj_ref->to_GFF3_format();
            
            print "$gff3_text\n";
            
        }
    }


    exit(0);
}

