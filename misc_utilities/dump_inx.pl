#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Gene_obj_indexer;

my $usage = "\n\nusage: $0 inx_file [key]\n\n"
    . "\tif no gene_id is given, all existing gene ids are provided via stdout\n\n";

my $inx_file = $ARGV[0] or die $usage;
my $key = $ARGV[1];

my $indexer = new TiedHash( { "use" => $inx_file } );

if (defined $key) {
    my $value = $indexer->get_value($key);
    print "[$key]\t$value\n";
}
else {
    ## provide the index listing:
    my @keys = $indexer->get_keys();
    foreach my $key (sort @keys) {
		my $value = $indexer->get_value($key);
		print "[$key]\t$value\n";
    }
}

exit(0);

