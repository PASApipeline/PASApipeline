#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Fasta_reader;
use GFF3_utils;
use Carp;


my $usage = "\n\nusage: $0 gff3_file\n\n";

my $gff3_file = $ARGV[0] or die $usage;

main: {

    my @all_gene_objs;

    my %had_UTR_trimmed;
    
    my $gene_obj_indexer_href = {};
    
    ## associate gene identifiers with contig id's.
    my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);
    
    foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
        
        print STDERR "// processing contig: $asmbl_id\n";
        
        my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
        
        my @gene_objs;
        
        foreach my $gene_id (@gene_ids) {
            my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
            push (@gene_objs, $gene_obj_ref);
        }
        
        push (@all_gene_objs, @gene_objs);

        #print STDERR "- got: " . scalar(@gene_objs) . " genes.\n";
        
        ## examine overlap
        for (my $i = 0; $i < $#gene_objs; $i++) {
            
            my $gene_i = $gene_objs[$i];
            
            for (my $j = $i + 1; $j <= $#gene_objs; $j++) {
                
                my $gene_j = $gene_objs[$j];
                
                
                if (&utr_overlaps_CDS($gene_i, $gene_j)) {
                    
                    &trim_utrs_remove_identical_CDS_isoforms($gene_i);
                    $had_UTR_trimmed{ $gene_i->{TU_feat_name} }++;
                    
                }
                
                if (&utr_overlaps_CDS($gene_j, $gene_i)) {
                    
                    &trim_utrs_remove_identical_CDS_isoforms($gene_j);
                    $had_UTR_trimmed{ $gene_j->{TU_feat_name} }++;
                    
                }
            }
        }
        
    }
    
    
    if (%had_UTR_trimmed) {
        
        #use Data::Dumper;
        #print Dumper(\%had_UTR_trimmed);
    
        foreach my $gene_obj (@all_gene_objs) {
            print $gene_obj->to_GFF3_format() . "\n";
        }
                
    }
    else {
        print STDERR "No UTRs required trimming.\n";
    }
    
    exit(0);
    
}


####
sub trim_utrs_remove_identical_CDS_isoforms {
    my ($gene) = @_;

    print STDERR "// removing UTRs from gene: " . $gene->toString() . "\n";

    &trim_utrs($gene);

    my %CDS_seen;
    my $redundancy_detected = 0;

    my @unique_isoforms;
    foreach my $g ($gene, $gene->get_additional_isoforms()) {
        
        my $cds_token = &get_CDS_token($g);
        if (exists ($CDS_seen{$cds_token})) {
            $redundancy_detected = 1;
        }
        else {
            $CDS_seen{$cds_token} = 1;
            push (@unique_isoforms, $g);
        }
    }

    if ($redundancy_detected) {
        print STDERR "REDUNDANCY DETECTED: " . $gene->toString() . "\n";
        
        my $parent_gene = shift @unique_isoforms;
        $parent_gene->delete_isoforms();
        if (@unique_isoforms) {
            $parent_gene->add_isoform(@unique_isoforms);
        }

        $parent_gene->refine_gene_object();
        
        print STDERR "*** REDUNDANCY FIXED: " . $parent_gene->toString();
        
        return($parent_gene);
    }
    else {
        return($gene);
    }

    
}


####
sub trim_utrs {
    my ($gene) = @_;

    foreach my $gene_obj ($gene, $gene->get_additional_isoforms()) {
        $gene_obj->trim_UTRs();
    }

    return;

}

    
####
sub utr_overlaps_CDS {
    my ($gene_A, $gene_B) = @_;

    #print "Checking overlap between: " . $gene_A->{TU_feat_name} . " and " . $gene_B->{TU_feat_name} . "\n";

    my ($A_lend, $A_rend) = sort {$a<=>$b} $gene_A->get_coords();
    my ($B_lend, $B_rend) = sort {$a<=>$b} $gene_B->get_coords();

    unless ($A_lend <= $B_rend && $A_rend >= $B_lend) {
        ## no overlap among genes
        return(0);
    }

    #print "Found overlap between: " . $gene_A->{TU_feat_name} . " and " . $gene_B->{TU_feat_name} . "\n";
    

    my @gene_A_utr_coords = &get_all_UTR_coordinates($gene_A);

    my @gene_B_CDS_coords = &get_all_CDS_coordinates($gene_B);

    foreach my $utr_coordset (@gene_A_utr_coords) {
        my ($utr_lend, $utr_rend) = sort {$a<=>$b} @$utr_coordset;
        
        foreach my $cds_coordset (@gene_B_CDS_coords) {
            my ($cds_lend, $cds_rend) = @$cds_coordset;
            
            if ($cds_lend <= $utr_rend && $cds_rend >= $utr_lend) {
                return(1);  # found overlap
            }
        }
    }

    return(0); # no overlap
    
}

####
sub get_all_UTR_coordinates {
    my ($gene_obj) = @_;

    my @utr_coords;

    foreach my $gene ($gene_obj, $gene_obj->get_additional_isoforms()) {

        my @p5_utr_coords = $gene->get_5prime_UTR_coords();
        
        if (@p5_utr_coords) {
            push (@utr_coords, @p5_utr_coords);
        }
        
        my @p3_utr_coords = $gene->get_3prime_UTR_coords();
        
        if (@p3_utr_coords) {
            push (@utr_coords, @p3_utr_coords);
        }
     
    }

    return(@utr_coords);
    
}

####
sub get_all_CDS_coordinates {
    my ($gene_obj) = @_;

    my @cds_coords;

    foreach my $gene ($gene_obj, $gene_obj->get_additional_isoforms()) {
        my @exons = $gene->get_exons();
        
        
        foreach my $exon (@exons) {
            if (my $cds = $exon->get_CDS_obj()) {
                my ($cds_end5, $cds_end3) = $cds->get_coords();
                push (@cds_coords, [$cds_end5, $cds_end3]);
            }
        }
        
    }
    return(@cds_coords);
}


####
sub get_CDS_token {
    my ($g) = @_;

    my $token = "";
    
    foreach my $exon ($g->get_exons()) {

        if (my $cds = $exon->get_CDS_exon_obj()) {

            my ($lend, $rend) = $cds->get_coords();
            $token .= "$lend-$rend,";
        }
    }

    chop $token; # remove trailing comma

    return($token);
}


    
    
        


