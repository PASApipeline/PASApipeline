#!/usr/bin/env perl

use strict;


while (<STDIN>) {
    if (/>(\S+)/) {
	print "$1\n";
    }
}


exit(0);

