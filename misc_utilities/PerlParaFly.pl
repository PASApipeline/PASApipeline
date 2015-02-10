#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use threads;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib/");
use Thread_helper;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use List::Util qw(shuffle);

my $usage = <<__EOUSAGE__;

##########################################################
#
# Usage: PerlParaFly (opts)
#
# Required: 
#   -c <str>              :filename containing list of bash-style commands to execute.
#   --CPU <int>            :number_of_threads
#
# Optional:
#   --shuffle              :randomly shuffles the command order. 
#   --failed_cmds <str>    :filename to capture failed commands.  default("FailedCommands")
#   -v                    :simple progress monitoring.
#
##########################################################

__EOUSAGE__

    ;


my $commands_file;
my $CPU;
my $shuffle_flag;
my $failed_cmds_file = "FailedCommands";
my $VERBOSE;
my $help;

&GetOptions ('c=s' => \$commands_file,
             'CPU=i' => \$CPU,
             'shuffle' => \$shuffle_flag,
             'failed_cmds=s' => \$failed_cmds_file,
             'v' => \$VERBOSE,
             'h' => \$help,
    );

if ($help) {
    die $usage;
}

unless ($CPU && $commands_file) {
    die $usage;
}

main: {

    my @cmds;
    open (my $fh, $commands_file) or die "Error, cannot open file $commands_file";
    while (<$fh>) {
        chomp;
        if (/\w/) {
            push (@cmds, $_);
        }
    }
    close $fh;
    
    if ($shuffle_flag) {
        @cmds = shuffle(@cmds);
    }


    my $thread_helper = new Thread_helper($CPU);
    
    my %thread_id_to_cmd;
    foreach my $cmd (@cmds) {
        
        my $thread = threads->create('process_cmd', $cmd);
        my $tid = $thread->tid();
        $thread_id_to_cmd{$tid} = $cmd;
        
        $thread_helper->add_thread($thread);
    }
    
    $thread_helper->wait_for_all_threads_to_complete();

    my @failures = $thread_helper->get_failed_threads();
    if (@failures) {
        
        print STDERR "Error, " . scalar(@failures) . " commands failed.  Writing to $failed_cmds_file\n";
        
        open (my $ofh, ">$failed_cmds_file") or die $!;
        foreach my $thread (@failures) {
            my $tid = $thread->tid();
            my $cmd = $thread_id_to_cmd{$tid};
            print $ofh "$cmd\n";
        }
        close $ofh;
        
    }
    else {
        print STDERR "Done. All commands completed successfully.\n";
    }
    
    exit(0);
}

####
sub process_cmd {
    my ($cmd) = @_;

    if ($VERBOSE) {
        print STDERR "CMD: $cmd\n";
    }

    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }

    return;
}

    
