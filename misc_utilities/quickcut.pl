#!/usr/bin/env perl 

use strict;
use warnings;
use lib ($ENV{EUK_MODULES});
use Nuc_translator;


my ($pos1, $pos2, @rest) = @ARGV;
my $orientation;
if ($pos1 < $pos2) {
    $orientation = '+';
} else {
    ($pos1, $pos2) = ($pos2, $pos1);
    $orientation = '-';
}

my ($header, $sequence);
while (<STDIN>) {
    chomp;
    if (/^>/) {
	$header = $_;
	next;
    }
    $_ =~ s/\s+//g;
    $sequence .= $_;
}

my $pull_length = ($pos2 - $pos1 + 1);
my $sub_sequence = substr ($sequence, ($pos1 - 1), $pull_length);
if ($orientation eq '-') {
    $sub_sequence = &reverse_complement($sub_sequence);
}
$sub_sequence =~ s/(\S{60})/$1\n/g;
chomp $sub_sequence;
print "$header $pos1 to $pos2 orient($orientation)\n$sub_sequence\n";


exit(0);
