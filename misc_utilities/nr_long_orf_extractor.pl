#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Gene_obj_indexer;
use CdbTools;

## reports only those non-redundant ORFs.

my $COMPLETE_ORF_ONLY_FLAG = 1;

my $usage = "\n\nusage: $0 gene_index_file min_ORF_length [GTF GenomeFastaDB]\n\n";

my $gene_inx = $ARGV[0] or die $usage;
my $min_orf_length = $ARGV[1] or die $usage;
my $format = ($ARGV[2] && $ARGV[2] eq 'GTF') ? "GTF" : "GFF3";
my $genome_fasta_db = $ARGV[3];
if ($format eq 'to_GTF_format' && ! -s $genome_fasta_db) {
    die "Error, need genome fasta file for GTF format.\n";
}

my $gene_obj_indexer = new Gene_obj_indexer( { "use" => $gene_inx } );

my %asmbl_id_to_geneobjs;

my @gene_objs;
my @gene_ids = $gene_obj_indexer->get_keys();


foreach my $gene_id (sort @gene_ids) {
    my $gene_obj_ref = $gene_obj_indexer->get_gene($gene_id);
    
    my $asmbl_id = $gene_obj_ref->{asmbl_id};
    
    die "Error, no asmbl_id " if ! $asmbl_id;
    
    my $gene_list_aref = $asmbl_id_to_geneobjs{$asmbl_id};
    unless (ref $gene_list_aref) {
        $gene_list_aref = $asmbl_id_to_geneobjs{$asmbl_id} = [];
    }
    push (@$gene_list_aref, $gene_obj_ref);
}

foreach my $asmbl_id (keys %asmbl_id_to_geneobjs) {    

    my $genome_seq;
    if ($format eq 'GTF') {
        $genome_seq = &cdbyank_linear($asmbl_id, $genome_fasta_db);
    }
    
    my %seen;
    my @gene_to_load;
    foreach my $gene (@{$asmbl_id_to_geneobjs{$asmbl_id}}) {
        
        my $cds_length = $gene->get_CDS_length();
        
        if ($cds_length / 3 < $min_orf_length) { next; }
        
        my @cds_coords;
        foreach my $exon ($gene->get_exons()) {
            my $cds = $exon->get_CDS_exon_obj();
            if ($cds) {
                push (@cds_coords, $cds->get_coords());
            }
        }
        
        unless (@cds_coords) { die "Error, no cds for " . $gene->toString(); }
        
        @cds_coords = sort {$a<=>$b} @cds_coords;
        my $coord_line = join (",", @cds_coords);
        if (! $seen{$coord_line}) {
            push (@gene_to_load, $gene);
            $seen{$coord_line} = 1;
        }
        
        if ($format eq 'GTF') {
           
            
            foreach my $isoform ($gene, $gene->get_additional_isoforms()) {
                
                $isoform->delete_isoforms(); # want GTF dumper to report only the isoform structure in unbundled form here.

                $isoform->create_all_sequence_types(\$genome_seq);
                
                my $gene_id = $isoform->{TU_feat_name};
                my $model_id = $isoform->{Model_feat_name};
                my $protein = $isoform->get_protein_sequence();
                my $cds_seq = $isoform->get_CDS_sequence();
        

                if ($COMPLETE_ORF_ONLY_FLAG) {
                    unless ($protein =~ /^M/ && $protein =~ /\*$/) { # got start and stop
                        print STDERR "ORF for $gene_id, $model_id is not complete. Skipping.\n";
                        next; 
                    }
                }
                
                print $isoform->to_GTF_format(\$genome_seq) . "\n";
     
                print "# $gene_id $model_id PROTEIN: $protein\n\n";
                print "# $gene_id $model_id CDS: $cds_seq\n\n";
            }
        }
        else {
            
            print $gene->to_GFF3_format() . "\n";
        }
        
    }
        
}


exit(0);
