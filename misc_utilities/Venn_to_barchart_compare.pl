#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 venn.txt\n\n";

my $venn = $ARGV[0] or die $usage;

my %class_to_att_count;

my %atts;

open (my $fh, $venn) or die "Error, cannot open file $venn";
while (<$fh>) {
	chomp;
	my ($feature, $class_list) = split(/\t/);

	my @classes = split(/,/, $class_list);
	
	my $att;
	if (scalar (@classes) == 1) {
		$att = "unique";
	}
	else {
		$att = "shared," . scalar(@classes);
	}
	
	$atts{$att}++;

	foreach my $class (@classes) {
		$class_to_att_count{$class}->{$att}++;
	}
}

close $fh;

my @sorted_atts = reverse sort keys %atts;

print "#class\t" . join("\t", @sorted_atts) . "\n";

foreach my $class (sort keys %class_to_att_count) {
	
	print $class;
	
	foreach my $att (@sorted_atts) {
		
		my $count = $class_to_att_count{$class}->{$att};
		print "\t$count";
	}
	
	print "\n";
}

exit(0);



		
	
	
