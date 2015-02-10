#!/usr/bin/env perl

## Pasa initialization

# setting PASAHOME env variable
# PASAHOME should either be the current working directory or the parent dir.

use FindBin;


BEGIN {
    
    ## Ensure that PASAHOME env var is set.
    ## Find PASAHOME by looking for the pasa_conf directory.
    
    unless ($ENV{PASAHOME}) {
        my $path = $FindBin::Bin;
        if (-d "$path/pasa_conf") {
            $ENV{PASAHOME} = $path;
        } elsif (-d "$path/../pasa_conf") {
            $ENV{PASAHOME} = "$path/../";
        } else {
            ## Search the Perl Lib Path.  
            foreach my $lib (@INC) {
                if (-d "$lib/pasa_conf") {
                    $ENV{PASAHOME} = $lib;
                    last;
                }
                if (-d "$lib/../pasa_conf") {
                    $ENV{PASAHOME} = "$lib/../";
                    last;
                }
            }
            
            unless ($ENV{PASAHOME}) {
                die "ERROR, cannot find the 'pasa_conf' directory. It should be the current directory or the parent dir.\n\n";
            }
        }
    }
}


use lib ("$ENV{PASAHOME}/PerlLib/", $ENV{PASAHOME});
use Pasa_conf;

1;
