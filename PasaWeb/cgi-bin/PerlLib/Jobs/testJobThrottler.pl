#!/usr/bin/env perl

use strict;
use warnings;

use Jobs::SysCmdJob;
use Jobs::Throttler;

my @jobs;
for (1..100) {
	push (@jobs, new Jobs::SysCmdJob("sleep ". int(rand(10))));
}

my $throttler = new Jobs::Throttler(4);
$throttler->launch(@jobs);

if (my @failed_jobs = $throttler->get_failed_jobs()) {
	die "Error, @failed_jobs failed\n";
}
else {
	print "All jobs succeeded.\n";
	exit(0);
}


