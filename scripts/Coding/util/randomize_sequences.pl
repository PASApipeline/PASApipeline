#!/usr/bin/env perl

use strict;
use warnings;

use List::Util qw (shuffle);
use FindBin;
use lib ("$FindBin::Bin/../../../PerlLib");
use Fasta_reader;

my $usage = "usage: $0 sequences.fa\n\n";

my $file = $ARGV[0] or die $usage;

my $seq_string = "";

my $fasta_reader = new Fasta_reader($file);
while (my $seq_obj = $fasta_reader->next()) {

	my $sequence = $seq_obj->get_sequence();
	$seq_string .= $sequence;
}

my @chars = split(//, $seq_string);

$seq_string = join("", shuffle(@chars));

$seq_string =~ s/(\S{60})/$1\n/g;

print ">randomized\n$seq_string\n";

exit(0);


