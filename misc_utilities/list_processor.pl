#!/usr/bin/env perl

use strict;

require "getopts.pl";

use vars qw ($opt_f $opt_S $opt_h $opt_n $opt_L $list $opt_d $DEBUG $splink @ids @cmds $cmd @x $asmbl);

$|=1;

&Getopts('dD:L:Sp:f:hn:');



if ($opt_h) {
    die <<_EOH_;

########################################################
#
# -L filename of list of asmbls to process
#
# -f filename of commands to run on each asmbl (place ASMBL where \$asmbl should be)
#
# -n number of processes to run simultaneously.
#
# -d debug
# -S verbose
# -h print this and exit.
#
#########################################################

_EOH_

}



if (length($opt_L) > 0) { $list = $opt_L; } else { &option_die;}
if (length($opt_d) > 0) { $DEBUG = 1; } else { $DEBUG = 0; }
my ($SEE) = ($opt_S) ? 1:0;
my ($number_processes_allowed) = ($opt_n) ? $opt_n : 1;
my $commands;
unless ( $commands = $opt_f) {die;}

# read in list of assemblies
open (LIST, $list) || die "Cant open $list";
while ($splink = <LIST>) {
    $splink =~ s/\s//g;
    print "SPLINK: $splink \n" if ($DEBUG);
    print "A: $splink\n" if ($DEBUG);
    push(@ids, $splink);
}

close(LIST);


## Read in command lines
print "Reading Commands\n" if $SEE;
open (COM, "$commands");
my $line;
while ($line = <COM>) {
    if ($line =~ /^\#/) {next;} #ignore commented out command lines.
    if ($line =~ /\w/) {
	print "Read Command: $line" if $SEE;
	chomp $line;
	push (@cmds, $line);
}
}
close COM;

my $counter = 0;

print "Total number of processes allowed at one time: $number_processes_allowed\n" if $SEE;

foreach $asmbl (@ids) {
    if ($asmbl) {
        print "\n------- Processing asmbl: $asmbl\n" if $SEE;
        my @copy_cmds = @cmds;
        foreach my $cmd (@copy_cmds) {
            $cmd =~ s/ASMBL/$asmbl/g;
            $cmd = "sleep 10" if $DEBUG;
            $counter++;
            if ($counter > $number_processes_allowed) {
                print "\tWaiting for a child to exit\n\n" if $SEE;
                wait();
            }
            print "CMD $cmd\n";
            my $whoami = fork();
            if (!$whoami) {
                #system "sleep 3"; #pause between launches for humanity's sake.
                print "\nRunning process $counter under child\n\n" if $SEE;
                system($cmd) if (!($DEBUG));
                print "Process $counter finished\n" if $SEE;
                exit();
            }
        }
    }
}
print "End of Parent Process\nWaiting for children to exit\n" if $SEE;
while (wait () > 0) {}


sub option_die {
    
    die;
}


