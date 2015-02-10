#!/usr/bin/env perl

use strict;
use warnings;
use lib ($ENV{EUK_MODULES});
use Gene_obj;
use GFF3_utils;
use Fasta_reader;
use Exons_to_geneobj;
use Data::Dumper;

my $usage = "usage: $0 old.gff3 new.gff3  genome.fasta\n\n";

my $old_gff3_file = $ARGV[0] or die $usage;
my $new_gff3_file = $ARGV[1] or die $usage;
my $genome_file = $ARGV[2] or die $usage;

my %STATS_COUNTER;

main: {
    
    ## parse old gff3 file
    my $old_gene_obj_indexer_href = {};
    my $old_contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($old_gff3_file, $old_gene_obj_indexer_href);
    
    ## parse new gff3 file
    my $new_gene_obj_indexer_href = {};
    my $new_contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($new_gff3_file, $new_gene_obj_indexer_href);
    
    my $fasta_reader = new Fasta_reader($genome_file);
    my %genome = $fasta_reader->retrieve_all_seqs_hash();
    
    
    foreach my $asmbl_id (sort keys %$old_contig_to_gene_list_href) {
        
        my $genome_seq = $genome{$asmbl_id} or die "Error, no genome sequence for scaffold $asmbl_id";
        
        my @gene_ids = @{$old_contig_to_gene_list_href->{$asmbl_id}};
        
        
        foreach my $gene_id (@gene_ids) {
            
            my $old_gene_obj = $old_gene_obj_indexer_href->{$gene_id};
            my $new_gene_obj = $new_gene_obj_indexer_href->{$gene_id};
            
            unless ($old_gene_obj && $new_gene_obj) {
                print STDERR "Warning, don't have both new and old versions of $gene_id\n";
                next;
            }

            $old_gene_obj->create_all_sequence_types(\$genome_seq);
            $new_gene_obj->create_all_sequence_types(\$genome_seq);

            my $old_protein = $old_gene_obj->get_protein_sequence();
            my $new_protein = $new_gene_obj->get_protein_sequence();
            
            if (length($old_protein) > length($new_protein)) {
                
                # print "Updated protein is truncated for $gene_id\n";
                
                &examine_for_retained_intron($old_gene_obj, $new_gene_obj, \$genome_seq);
                
            }
            
            
            
        }
    }

    
    my $percent_fixed = $STATS_COUNTER{FIXED}/$STATS_COUNTER{RI} * 100;

    print STDERR "Fixed " . sprintf("%.2f", $percent_fixed) . "% (" . $STATS_COUNTER{FIXED} . "/" . $STATS_COUNTER{RI} . ") retained-intron containing genes.\n";
    
        
    exit(0);
    

}

####
sub examine_for_retained_intron {
    my ($old_gene_obj, $new_gene_obj, $genome_seq_sref) = @_;
    
    my $gene_id = $old_gene_obj->{TU_feat_name};

    my @old_intron_coords = $old_gene_obj->get_intron_coordinates();

    unless (@old_intron_coords) {
        ## no introns
        return;
    }

    my $found_retained_intron_flag = 0;

    my @revisions;

    foreach my $exon ($new_gene_obj->get_exons()) {

        my ($exon_lend, $exon_rend) = sort {$a<=>$b} $exon->get_coords();
        
        my $struct = { exon_coords => [$exon_lend, $exon_rend],
                       retained_introns => [],
                   };
        
        push (@revisions, $struct);
        
        
        foreach my $intron_coordset (@old_intron_coords) {

            my ($intron_lend, $intron_rend) = sort {$a<=>$b} @$intron_coordset;
            
            if ($exon_lend < $intron_lend && $intron_rend < $exon_rend) {
                
                # print join("\t", $gene_id, "exon $exon_lend-$exon_rend harbors", "intron $intron_lend-$intron_rend") . "\n";
                
                $found_retained_intron_flag = 1;
                
                push (@{$struct->{retained_introns}}, [$intron_lend, $intron_rend]);
            
            }
        }
    }

    if ($found_retained_intron_flag) {
        
        $STATS_COUNTER{RI}++;
        
        
        my $altered_gene_obj = &reinsert_spliced_introns($new_gene_obj, \@revisions, $genome_seq_sref);
        
        $altered_gene_obj->create_all_sequence_types($genome_seq_sref);
        
        my $altered_prot_seq = $altered_gene_obj->get_protein_sequence();
        
        my $old_prot_seq = $old_gene_obj->get_protein_sequence();
        
        if (length($altered_prot_seq) >= length($old_prot_seq) ) {

            
            $altered_gene_obj->{asmbl_id} = $old_gene_obj->{asmbl_id};
            $altered_gene_obj->{TU_feat_name} = $old_gene_obj->{TU_feat_name};
            $altered_gene_obj->{Model_feat_name} = $old_gene_obj->{Model_feat_name};
            $altered_gene_obj->{com_name} = $old_gene_obj->{com_name};
            
            # print "Restoring introns lengthened the protein sequence.\n";
            print $altered_gene_obj->to_GFF3_format() . "\n#$altered_prot_seq\n\n";
            
            $STATS_COUNTER{FIXED}++;
        }
    }
    

    return;
}


####
sub reinsert_spliced_introns {
    my ($new_gene_obj, $revisions_aref, $genome_seq_sref) = @_;

    my $orient = $new_gene_obj->get_orientation();
        
    my @exon_coordsets = ();

    foreach my $revision_struct (@$revisions_aref) {

        my $exon_coordset_aref = $revision_struct->{exon_coords};
        my $retained_introns_aref = $revision_struct->{retained_introns};

        if (@$retained_introns_aref) {

            ## split the exon up.
            my @intron_coords = sort {$a->[0]<=>$b->[0]} @$retained_introns_aref;
            
            my $exon_lend = $exon_coordset_aref->[0];

            foreach my $intron (@intron_coords) {
                my ($intron_lend, $intron_rend) = @$intron;
                
                my $new_exon_rend = $intron_lend - 1;
                push (@exon_coordsets, [$exon_lend, $new_exon_rend]);
                
                $exon_lend = $intron_rend + 1;
            }
            push (@exon_coordsets, [$exon_lend, $exon_coordset_aref->[1]]);
        }
        else {
            ## keep original exon
            push (@exon_coordsets, $exon_coordset_aref);
        }

    }

    my %new_exon_coords;
    foreach my $coordset (@exon_coordsets) {
        
        my ($lend, $rend) = @$coordset;
       
        my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);
        
        $new_exon_coords{$end5} = $end3;
        
    }
    
    #print "New exon coords: " . Dumper(\%new_exon_coords);
    
    my $altered_gene_obj = &Exons_to_geneobj::create_gene_obj(\%new_exon_coords, $genome_seq_sref);

    #print "Altered gene: " . $altered_gene_obj->toString();
    

    return ($altered_gene_obj);
}


                
