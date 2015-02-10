#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\nusage: $0 compreh_init_build.details RSEM.isoforms.results\n\n";

my $query_to_OS_summary_file = $ARGV[0] or die $usage;
my $RSEM_isoforms_file = $ARGV[1] or die $usage;

main: {

    my %trans_to_cat;
    { 
        open (my $fh, $query_to_OS_summary_file) or die $!;
        while (<$fh>) {
            chomp;
            my ($gene, $trans, $class) = split(/\t/);
            
            $trans_to_cat{$trans} = $class;
            
        }
        close $fh;
    }
    
    my %class_to_expr_sum;
    my $total_fpkm = 0;
    {
        
        open (my $fh, $RSEM_isoforms_file) or die $!;
        my $header = <$fh>;
        while (<$fh>) {
            chomp;
            my @x = split(/\t/);
            my $acc = $x[0];
            my $fpkm = $x[6];
            
            my $class = $trans_to_cat{$acc} or die "Error, no class assigned to acc: $acc ";

            $class_to_expr_sum{$class} += $fpkm;
            
            $total_fpkm += $fpkm;
            
        }
        close $fh;
        

    }

    foreach my $class (sort keys %class_to_expr_sum) {
        
        my $count = $class_to_expr_sum{$class};
        my $pct = sprintf("%.2f", $count/$total_fpkm*100);

        print join("\t", $class, $count, $pct) . "\n";
    }

    exit(0);
}


    
    
