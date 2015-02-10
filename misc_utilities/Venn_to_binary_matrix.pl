#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 venn.txt\n\n";

my $venn = $ARGV[0] or die $usage;

=R_code_show_dendrogram

 data = read.table("venn.txt.bmatrix", com='', row.names=1, header=T)
 d = dist(as.matrix(t(data)), method='binary')
 hc = hclust(d)
 plot(hc)


=cut



my %atts;
my %entries;

open (my $fh, $venn) or die $!;
while (<$fh>) {
    chomp;
    my ($entry, $att_types) = split(/\t/);
    foreach my $att (split(/,/, $att_types)) {

        $atts{$att} = 1;
        $entries{$entry}->{$att} = 1;
    }

}


my @atts = sort keys %atts;
print "#entry\t" . join("\t", @atts) . "\n";
foreach my $entry (keys %entries) {
    print $entry;
    foreach my $att (@atts) {
        my $val = $entries{$entry}->{$att} || 0;
        print "\t$val";
    }
    print "\n";
}

exit(0);


    
