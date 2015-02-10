#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 cmds.list\n\n";

my $cmds_file = $ARGV[0] or die $usage;


main: {
	
	open (my $fh, $cmds_file) or die "Error, cannot open file $cmds_file";
	while (<$fh>) {
		chomp;
		
		if (/^CMD: /) {
			s/^CMD: //;
			&process_cmd($_);
		}
	}

	exit(0);

}

####
sub process_cmd {
	my ($cmd) = @_;

	print "EXECUTING: $cmd\n";

	my $ret = system($cmd);

	if ($ret) {
		die "Error, cmd: $cmd died with ret $ret";
	}

	return;
}


