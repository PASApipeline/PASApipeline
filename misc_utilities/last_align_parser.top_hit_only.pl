#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\nusage: $0 lastal.parsed\n\n";


my $lastal_parsed_file = $ARGV[0] or die $usage;

main: {


    my %data;

    open (my $fh, $lastal_parsed_file) or die $!;
    while (<$fh>) {
        my $line = $_;
        if (/^\#/) { next; }
        chomp;
        my @x = split(/\t/);
        
        my $read_acc = $x[0];
        my $coords = $x[3];
        my $match_len = 0;
        my $num_aligned_chars = $x[4];
        my $per_id = $x[5];

        my $match_score = $per_id * $num_aligned_chars;


        if ( (! exists $data{$read_acc}) 
             ||
             $match_score > $data{$read_acc}->{score}) {

            
            $data{$read_acc} = { score => $match_score,
                                 line => $line,
                             };
        }

    }
    close $fh;

    foreach my $struct (values %data) {
        print $struct->{line};
    }

    exit(0);
}
