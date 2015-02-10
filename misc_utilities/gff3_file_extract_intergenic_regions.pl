#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use CdbTools;
use GFF3_utils;
use Carp;
use Nuc_translator;
use Longest_orf;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case bundling);


my $usage = <<__EOUSAGE__;

#############################################################################################################
#
#  --gff3 <string>      annotations in GFF3 format
#
#  --trim <int>         percentile of each end of the intergenic feature length distribution to exclude (default: 20%)
#
#  --flank <int>        length to add to the end of each gene (account for possible UTR lengths) (default: 0);
#
##############################################################################################################

__EOUSAGE__

    ;


my $gff3_file;
my $trim_percent = 20;
my $flank_add = 0;

&GetOptions ( "gff3=s" => \$gff3_file,
              "trim=i" => \$trim_percent,
              "flank=i" => \$flank_add,
              );


unless ($gff3_file) {
    die $usage;
}


my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

my @intergenic_regions;

foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
    my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
    my @genes;

    foreach my $gene_id (@gene_ids) {
        my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
        
        my ($lend, $rend) = sort {$a<=>$b} $gene_obj_ref->get_coords();

        $lend -= $flank_add;
        $rend += $flank_add;
        
        my $struct = { lend => $lend,
                       rend => $rend,
                   }; 
        
        push (@genes, $struct);
    }

    @genes = sort {$a->{lend}<=>$b->{lend}} @genes;

    my $prev_gene = shift @genes;
    my $prev_rend = $prev_gene->{rend};

    foreach my $gene (@genes) {
                
        my $curr_lend = $gene->{lend};
        my $curr_rend = $gene->{rend};

        if ($curr_lend > $prev_rend) {
            
            my $intergenic_lend = $prev_rend + 1;
            my $intergenic_rend = $curr_lend - 1;

            my $intergene = { lend => $intergenic_lend,
                              rend => $intergenic_rend,
                              length => $intergenic_rend - $intergenic_lend + 1,
                              scaffold => $asmbl_id,
                          };

            push (@intergenic_regions, $intergene);
        }
        
        if ($curr_rend > $prev_rend) {
            $prev_rend = $curr_rend;
        }
        
    }

}


@intergenic_regions = sort {$a->{length}<=>$b->{length}} @intergenic_regions;


if ($trim_percent > 0) {
    # trim it
    my $num_intergenic_regions = scalar(@intergenic_regions);
    
    my $index_low_trim = int ($num_intergenic_regions * ($trim_percent/100));
    my $index_high_trim = $num_intergenic_regions - $index_low_trim;
    
    @intergenic_regions = @intergenic_regions[$index_low_trim..$index_high_trim];
}


## output intergenic regions in GTF format

foreach my $intergenic_region (@intergenic_regions) {
    
    my $scaffold = $intergenic_region->{scaffold};
    my $lend = $intergenic_region->{lend};
    my $rend = $intergenic_region->{rend};

    my $gene_obj = new Gene_obj();
    my $coords_href = { $lend => $rend };
    
    $gene_obj->populate_gene_object({}, $coords_href);

    $gene_obj->{com_name} = "intergene.$scaffold.$lend-$rend";
    
    $gene_obj->{asmbl_id} = $scaffold;
    $gene_obj->{TU_feat_name} = "intergene.$scaffold.$lend-$rend";
    $gene_obj->{Model_feat_name} = "m.intergene.$scaffold.$lend-$rend";
    

    print $gene_obj->to_transcript_GTF_format() . "\n";

    

}


exit(0);


        
                       
