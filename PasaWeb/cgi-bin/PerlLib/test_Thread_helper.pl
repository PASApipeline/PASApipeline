#!/usr/bin/env perl

use FindBin;
use lib ("$FindBin::Bin");

use threads;
use Thread_helper;
use Data::Dumper;

my $num_simultaneous_threads = 10;
my $thread_helper = new Thread_helper($num_simultaneous_threads);

for my $num (1..1000) {

    print STDERR "$num\n";

    $thread_helper->wait_for_open_thread();

    my $thread = threads->create(process_cmd, 'sleep', int(rand(10)));
    $thread_helper->add_thread($thread);
}
$thread_helper->wait_for_all_threads_to_complete();

my @failures = $thread_helper->get_failed_threads();
if (@failures) {
        ## examine them...   these are the same threads created above, use use the threads api to access info about them
        ## such as error messages
    print STDERR Dumper(\@failures) . " failed.";
    exit(1);
}
else {
    ## all good!
    
    print STDERR "Done.\n";
    exit(0);
}


####
sub process_cmd {
    my (@params) = @_;


    print STDERR "CMD: @params\n";
    my $ret = system(@params);

    if ($ret) {
        die "Error, CMD: @params died with ret $ret";
    }

    return;
}
