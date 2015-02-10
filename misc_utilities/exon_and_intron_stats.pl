#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/PerlLib");
use Gene_obj;
use Gene_obj_indexer;

my $usage = "\n\nusage: $0 gene_index_file [gene_id]\n\n"
    . "\tif no gene_id is given, all existing gene ids are provided via stdout\n\n";

my $gene_inx = $ARGV[0] or die $usage;
my $gene_id = $ARGV[1];


my $gene_obj_indexer = new Gene_obj_indexer( { "use" => $gene_inx } );

my $exon_count = 0;
my $intron_count = 0;
my $exon_lengths_sum = 0;
my $intron_lengths_sum = 0;

my $single_exon_gene_count = 0;
my $multi_exon_gene_count = 0;

if ($gene_id) {
    my $gene_obj = $gene_obj_indexer->get_gene($gene_id);
    &add_gene_to_stats($gene_obj);
}
else {
    ## provide the index listing:
    my @gene_ids = $gene_obj_indexer->get_keys();
    foreach my $gene_id (sort @gene_ids) {
        my $gene_obj = $gene_obj_indexer->get_gene($gene_id);
        &add_gene_to_stats($gene_obj);
    }
}


my $average_cds_exon_length = sprintf ("%.2f", $exon_lengths_sum / $exon_count);
my $average_intron_length = sprintf ("%.2f", $intron_lengths_sum / $intron_count);

my $total_num_genes = $single_exon_gene_count + $multi_exon_gene_count;
my $percent_single_exon = $single_exon_gene_count / $total_num_genes;

printf "$single_exon_gene_count single exon genes = %.2f %\n", $percent_single_exon * 100;
print "$multi_exon_gene_count multi_exon_genes\n";

print "CDS-exons:  count: $exon_count  average length: $average_cds_exon_length\n";
print "Introns:   count: $intron_count average length: $average_intron_length\n";



exit(0);


####
sub add_gene_to_stats {
    my ($gene_obj) = @_;
    
    my @exons = $gene_obj->get_exons();
    foreach my $exon (@exons) {
        if (my $cds = $exon->get_CDS_obj()) {
            my $length = $cds->length();
            $exon_count++;
            $exon_lengths_sum += $length;
        }
    }

    my @introns = $gene_obj->get_intron_coordinates();
    
    if (scalar @introns == 0) {
        $single_exon_gene_count++;
    }
    else {
        $multi_exon_gene_count++;
    }


    foreach my $intron (@introns) {
        my ($lend, $rend) = sort {$a<=>$b} @$intron;
        my $length = $rend - $lend + 1;
        $intron_count++;
        $intron_lengths_sum += $length;
    }
}


