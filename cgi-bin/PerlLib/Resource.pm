#!/usr/bin/env perl

package Resource;

use strict;
use warnings;

use FindBin;



## search the @INC paths for the resource and return the directory containing the file or package. (the @INC entry itself).
####
sub find_resource {
	my $filename_or_package = @_;


	$filename_or_package =~ s|::|/|g;

	foreach my $dir (@INC) {
		if (-e "$dir/$filename_or_package" || "$dir/$filename_or_package.pm") {
			return ($dir);
		}
	}

	return;
}


1;

