#!/usr/bin/env perl

use strict;
use warnings;

use ColorGradient;

my $num_colors = $ARGV[0] or die "usage: $0 num_colors\n\n";

my @colors = &ColorGradient::get_RGB_gradient($num_colors);

my @hex_colors = &ColorGradient::convert_RGB_hex(@colors);

print "<table>\n";
foreach my $color (@hex_colors) {
    print "<tr><td bgcolor=\'$color\'>$color</td></tr>\n";
	
}
print "</table>\n";

exit(0);

