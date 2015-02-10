#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin");
use XMFA_parser;


my $usage = "usage: $0 alignments.xmfa\n\n";

my $alignments_xmfa_file = $ARGV[0] or die $usage;


main: {
    
    
    my $xmfa_parser = new XMFA_parser($alignments_xmfa_file);

    $xmfa_parser->print();

    exit(0);
}

