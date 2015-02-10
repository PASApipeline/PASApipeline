#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Gene_obj_indexer;
use CdbTools;
use GFF3_utils;
use Carp;

$|++;

my $usage = "\n\nusage: $0 gff3_file [ gff3_file, ... ]\n\n"
    . "\tgenes are indexed by model feat_name.\n\n";


my @gff3_files = @ARGV;
unless (@gff3_files) { die $usage; }

my $index_file = "gene_structures.inx";
if (scalar @gff3_files == 1) {
    $index_file = $gff3_files[0] . ".inx";
}

my $gene_obj_indexer = new Gene_obj_indexer( { "create" => $index_file } );

my %seen; #track gene_ids already visited

foreach my $gff3_file (@gff3_files) {
    
    &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer);
    
}

exit(0);
