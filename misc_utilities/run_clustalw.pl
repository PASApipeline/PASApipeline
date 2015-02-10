#!/usr/local/bin/perl

use strict;

my $usage = "usage: $0 multi_fasta_file";

my $fasta_file = $ARGV[0] or die "\n$usage\n\n";

$fasta_file =~ s/^[^=]+=//;

my $cmd = "clustalw -output=gcg -infile=${fasta_file} -outfile=${fasta_file}.msf";
#my $cmd = "clustalw  -infile=$fasta_file -OUTORDER=INPUT > /dev/null ";

print "CMD: $cmd\n";
system ($cmd);

