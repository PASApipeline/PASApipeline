#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 analysis.txt blat_alignments.psl [token='|BLAT_mapped']\n\n";

my $analysis_txt = $ARGV[0] or die $usage;
my $blat_alignments = $ARGV[1] or die $usage;
my $token = $ARGV[2] || "|BLAT_mapped";

main: {

    my %blat_mapped;
    {
        open (my $fh, $blat_alignments) or die $!;
        while (<$fh>) {
            chomp;
            my @x = split(/\t/);
            my $acc = $x[9];
            $blat_mapped{$acc} = 1 if $acc;
        }
        close $fh;
    }

    open (my $fh, $analysis_txt) or die $!;
    while (<$fh>) {
        
        my @x = split(/\t/);
        if ($blat_mapped{$x[1]}) {
            $x[0] .= $token;
        }
        
        print join("\t", @x);
    }
    close $fh;

    exit(0);
}


        
