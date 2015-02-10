#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 keys1 file2 [field_no]\n\n";
my $fileA = $ARGV[0] or die $usage;
my $fileB = $ARGV[1] or die $usage;
my $field_no = $ARGV[2] || 0;

my %tokens;
open (my $fh, $fileA) or die "Error, cannot open $fileA ";
while (<$fh>) {
    s/^\s+//; #trim leading whitespace
    my @x = split (/\s+/);
    $tokens{$x[0]} = 1;
}
close $fh;

open ($fh, $fileB) or die "Error, cannot open $fileB ";
while (<$fh>) {
    my @x = split (/\s+/);
    if ($tokens{$x[$field_no]}) {
        print;
    }
}
close $fh;


exit(0);


        
   
