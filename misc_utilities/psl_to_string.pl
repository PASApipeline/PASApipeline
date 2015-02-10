#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use PSL_parser;

my $usage = "usage: $0 pslfile|STDIN\n\n";

my $file = $ARGV[0] || *STDIN{IO};

main: {
	
	my $psl_parser = new PSL_parser($file);
	while (my $pe = $psl_parser->get_next()) {
		print $pe->toString() . "\n";
	}
	
	exit(0);
}



	
