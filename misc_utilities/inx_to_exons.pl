#!/usr/local/bin/perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Gene_obj_indexer;
use GFF3_utils;
use Carp;

$|++;

my $usage = "\n\nusage: $0 file.inx\n\n";

my $inx_file = $ARGV[0] or die $usage;

unless (-s $inx_file) {
    die $usage;
}

my $gene_obj_indexer = new Gene_obj_indexer( { "use" => $inx_file } );

my @gene_ids = $gene_obj_indexer->get_keys();
foreach my $gene_id (@gene_ids) {
    
    

    my $gene_obj = $gene_obj_indexer->get_gene($gene_id);
    
    my $contig_id = $gene_obj->{asmbl_id};
    
    my @intron_coords = $gene_obj->get_intron_coordinates();
    my $strand = $gene_obj->get_orientation();
    my $model_id = $gene_obj->{Model_feat_name};
    my $source = $gene_obj->{source} || ".";
    
    my @exons = $gene_obj->get_exons();
    ## be sure they're sorted properly.
    @exons = sort {$a->{end5}<=>$b->{end5}} @exons;
    if ($strand eq '-') {
        @exons = reverse @exons;
    }


    if (scalar @exons > 1) {
        my $initial_exon =  shift @exons;
        $initial_exon->{type} = "Initial";
        my $last_exon = pop @exons;
        $last_exon->{type} = "Terminal";
        foreach my $remaining_exon (@exons) {
            $remaining_exon->{type} = "Internal";
        }

        unshift (@exons, $initial_exon);
        push (@exons, $last_exon);
    }
    else {
        $exons[0]->{type} = "Single";
    }
    
    foreach my $exon (@exons) {
        my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
        print "$contig_id\t$source\t$model_id\t[$strand]\t" . $exon->{type} . "\t$lend\t$rend\n";
    }
    
}

exit(0);



