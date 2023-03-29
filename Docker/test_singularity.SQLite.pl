#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use File::Spec;

my $base_dir = File::Spec->rel2abs("$FindBin::Bin/../");



my $version = `cat VERSION.txt`;
chomp $version;

my $sample_data_path = "$base_dir/sample_data";

chdir($sample_data_path) or die $!;

my $cmd = "singularity exec -e -B $base_dir  ../Docker/pasapipeline.v${version}.simg ./runMe.SQLite.sh";

print "CMD: $cmd\n";

my $ret = system($cmd);

exit ($ret);


