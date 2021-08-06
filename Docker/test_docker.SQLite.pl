#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

my $base_dir = "$FindBin::Bin/../";


my $cmd = "docker run --rm -it -v /tmp:/tmp -v $base_dir:$base_dir pasapipeline/pasapipeline:latest "
    . " bash -c 'cd /$base_dir/sample_data && ./runMe.SQLite.sh' ";

print "CMD: $cmd\n";

my $ret = system($cmd);

exit ($ret);


