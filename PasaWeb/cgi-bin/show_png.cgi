#!/usr/bin/env perl

use strict;

use CGI;
use CGI::Carp qw(fatalsToBrowser);

my $cgi = new CGI();
my $image = $cgi->param('image');
my $keep_image = $cgi->param('keep_image');

print $cgi->header(-type=>'image/gif');

open (FILE, $image) or die $!;
while (<FILE>) {
    print;
}
close FILE;

unless ($keep_image) {
    unlink ($image);
}

exit(0);
