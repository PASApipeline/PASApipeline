#!/usr/bin/env perl

use strict;
use warnings;

while (<STDIN>) {
    chomp;
    
    if (/error/i) { next;}
    
    my @x = split (/\t/);
    

    my ($query_acc, $query_len, $db_acc, $per_id, $db_len, $evalue) = ($x[0], $x[2], $x[5], $x[10], $x[18], $x[19]);
    
    unless ($per_id > 1) {
        next;
    }
    
    my ($match_len) = ($x[9] - $x[8]) + 1;
    
    my $percent_db_len = sprintf ("%.2f", $match_len / $db_len * 100);
    my $percent_query_len = sprintf ("%.2f", $match_len / $query_len * 100);

    print "$query_acc\t$db_acc\tper_id: $per_id\tper_query_len: $percent_query_len\tper_db_len: $percent_db_len\tevalue: $evalue\n";
}


