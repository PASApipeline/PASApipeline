#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../");
use Rlite;

## The McDonald and Kreitman test as desribed in:
##   Adaptive protein evolution at the Adh locus in Drosophila   Nature 1991


print <<_EOEXAMPLE_;

## the contingency table:
#                
#                            Fixed (between species)           Polymorphic (within species)
#         Replacement         7                                 2
#         Synonymous         17                                42
#
#       G=7.43, P=0.006     G-test of independence using the Williams correction for continuity
#

_EOEXAMPLE_

	;

my @contingency_table = ( [7, 2], [17, 42] );

my ($Gval, $prob) = &Rlite::G_test_matrix(\@contingency_table, "williams");

print "Gval: $Gval\tProb: $prob\n";

exit(0);

