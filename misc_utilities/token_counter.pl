#!/usr/bin/env perl

use strict;

my %tokens;

while (<STDIN>) {
    chomp;
    $tokens{$_}++;
}

foreach my $token (reverse sort {$tokens{$a}<=>$tokens{$b}} keys %tokens) {
    print $token . "\t$tokens{$token}\n";
}


