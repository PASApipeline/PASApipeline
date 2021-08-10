#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

my $base_dir = "$FindBin::Bin/../";


my $version = `cat VERSION.txt`;
chomp $version;

my $cmd = "docker run --rm -it -v $base_dir:$base_dir pasapipeline/pasapipeline:$version "
    . " bash -c 'service mysql start && cd $base_dir/sample_data && ./runMe.MySQL.sh' ";

print "CMD: $cmd\n";

my $ret = system($cmd);

exit ($ret);


