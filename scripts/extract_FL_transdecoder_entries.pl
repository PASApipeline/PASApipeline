#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 transdecoder.gff3\n\n";

my $td_gff3 = $ARGV[0] or die $usage;

main: {
    
    open (my $fh, $td_gff3) or die $!;
    while (<$fh>) {
        chomp;
        unless (/\w/) { next; }
        my @x = split(/\t/);
        if ($x[2] eq "gene" && $x[8] =~ /complete/) {
            print "$x[0]\n";
        }
    }
    close $fh;
    
    exit(0);
}

