#!/usr/bin/env perl

use strict;
use warnings;

use Carp;
use FindBin;
use lib ("$FindBin::Bin/");

use PSL_parser;

my $usage = "usage: $0 file.psl\n\n";

my $psl_file = $ARGV[0] or die $usage;

main: {

	my $psl_parser = new PSL_parser($psl_file);
	
	while (my $psl_entry = $psl_parser->get_next()) {

		print $psl_entry->toString() . "\n";

	}


	exit(0);
}


