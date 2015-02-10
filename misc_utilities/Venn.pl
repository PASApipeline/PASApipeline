#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;

my $usage = "usage: $0 file1 file2 file3 ...\n\n";

my @atts = @ARGV;
unless (@atts) {
    die $usage;
}

my %ele_to_att_list;

foreach my $file (@atts) {
    open (my $fh, $file) or die "Error, cannot open $file\n";
    
	$file = basename($file);
	
	while (<$fh>) {
        chomp;
        s/^\s+|\s+$//g;
        my $ele = $_;
        
        my $att_list_aref = $ele_to_att_list{$ele};
        unless (ref $att_list_aref) {
            $att_list_aref = $ele_to_att_list{$ele} = [];
        }
        push (@$att_list_aref, $file) unless (grep {$_ eq $file} @$att_list_aref); # don't add already existing attribute
                
    }
    close $fh;
}

foreach my $ele (sort keys %ele_to_att_list) {
    my $att_list_aref = $ele_to_att_list{$ele};
    print "$ele\t" . join (',', @$att_list_aref) . "\n";
}

exit(0);

