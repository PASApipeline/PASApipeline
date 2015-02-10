#!/usr/bin/env perl

use strict;
use warnings;
use lib ($ENV{EUK_MODULES});
use Gene_obj;


my $usage = "usage: $0 file.gff\n\n";

my $file = $ARGV[0] or die $usage;

main: {

    my %gene_coords;

    open (my $fh, $file) or die "Error, cannot open file $file";
    while (<$fh>) {
        
        chomp;
        my @x = split(/\t/);
        my $chr = $x[0];
        my $type = $x[2];
        my $lend = $x[3];
        my $rend = $x[4];
        my $orient = $x[6];
        my $info = $x[8];

        my $id;
        if ($info =~ /Parent=([^;\s]+)/) {
            $id = $1;
        }
        elsif ($info =~ /ID=([^;]+)/) {
            $id = $1;
        }
        
        unless ($id) {
            die "Error, no ID parsed for $_";
        }

        my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);

        $gene_coords{$id}->{$type}->{$end5} = $end3;
        $gene_coords{$id}->{chr} = $chr;
        
    }

    foreach my $gene_id (keys %gene_coords) {

        my $gene_obj = new Gene_obj();
        
        my $chr = $gene_coords{$gene_id}->{chr};
        
        my $cds_href = $gene_coords{$gene_id}->{CDS};

        $gene_obj->populate_gene_object($cds_href, $cds_href);
        
        $gene_obj->{asmbl_id} = $chr;
        $gene_obj->{TU_feat_name} = $gene_id;
        $gene_obj->{com_name} = $gene_id;
        $gene_obj->{Model_feat_name} = "m.$gene_id";
        

        print $gene_obj->to_GFF3_format() . "\n";
    }


    exit(0);
}

