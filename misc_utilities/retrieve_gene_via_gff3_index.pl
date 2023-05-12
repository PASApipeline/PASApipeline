#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Gene_obj_indexer;

my $usage = "\n\nusage: $0 gene_index_file [gene_id|file_with_gene_id_list]\n\n"
    . "\tif no gene_id is given, all existing gene ids are provided via stdout\n\n";

my $gene_inx = $ARGV[0] or die $usage;
my $gene_id = $ARGV[1];


my $gene_obj_indexer = new Gene_obj_indexer( { "use" => $gene_inx } );


if ($gene_id) {
    if (-f $gene_id) {
        # get each from file
        my $gene_id_file = $gene_id;
        open (my $fh, $gene_id_file) or die "Error, cannot open file $gene_id_file";
        while(<$fh>) {
            my $gene_id = $_;
            chomp $gene_id;
            my $gene_obj = $gene_obj_indexer->get_gene($gene_id);
            print $gene_obj->to_GFF3_format();
        }
    }
    else {
        my $gene_obj = $gene_obj_indexer->get_gene($gene_id);
        print $gene_obj->to_GFF3_format();
    }
}
else {
    ## provide the index listing:
    my @gene_ids = $gene_obj_indexer->get_keys();
    foreach my $gene_id (sort @gene_ids) {
        print "[$gene_id]\n";
    }
}

exit(0);

