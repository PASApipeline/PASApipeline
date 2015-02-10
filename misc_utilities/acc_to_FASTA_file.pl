#!/usr/bin/env perl

use strict;
use warnings;
use lib ($ENV{EUK_MODULES});
use CdbTools;

my $usage = "\nusage: $0 acc_listing fastaDB\n\n";

my $acc_listing = $ARGV[0] or die $usage;
my $fasta_db = $ARGV[1] or die $usage;

open (ACCS, "$acc_listing") or die "Can't open $acc_listing";

while (<ACCS>) {
    #print;
    chomp;
    my ($acc, $header) = split (/\s+/, $_, 2);
    
    unless (defined $header) {
        $header = "";
    }
    
	my $fasta_entry = &cdbyank($acc, $fasta_db);
	
	my @components = split (/\n/, $fasta_entry);
	my $fasta_header = shift @components;
	
	print "$fasta_header $header\n" . join ("\n", @components) . "\n";
	
}

exit(0);



