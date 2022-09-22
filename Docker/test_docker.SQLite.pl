#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

my $base_dir = "$FindBin::Bin/../";


my $version = `cat VERSION.txt`;
chomp $version;

if (! -d "/tmp") {
    mkdir "/tmp";
}

my $cmd = "docker run --rm -it -v `pwd`/tmp:/tmp -v $base_dir:$base_dir pasapipeline/pasapipeline:$version "
    . " bash -c 'cd $base_dir/sample_data && ./runMe.SQLite.sh' ";

print "CMD: $cmd\n";

my $ret = system($cmd);

exit ($ret);


