#!/usr/bin/env perl

use strict;
use warnings;

my %all_objects;
my %sets;

my $file = $ARGV[0] or die "usage: $0 Venn.output\n\n";

open (my $fh, $file) or die "Error, cannot open file $file";

## populate classes
while (<$fh>) {
	chomp;
	my ($acc, $venn_list) = split (/\s+/);
	my @atts = split (/,/, $venn_list);
	foreach my $att (@atts) {
		push (@{$sets{$att}}, $acc);
	}
	
	$all_objects{$acc} = 1;
}
close $fh;

print "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>

<ClassificationTree version=\"1.0\">

<ObjectSet>\n";

foreach my $acc (keys %all_objects) {
	print "<Object ID=\"$acc\"><Name>$acc</Name><Location>none</Location></Object>\n";
}
print "</ObjectSet>

<ClassificationSet>\n";

print "<Classification ID=\"root\">\n"
	. "<Name>16Sregions</Name>\n"
	. "<Objects objectIDs=\"" . join (" ", keys %all_objects) . "\"/>\n"
	. "</Classification>\n\n";

foreach my $class (keys %sets) {
	print "<Classification ID=\"$class\">\n";
	print "<Name>$class</Name>\n";
	print "<SuperClass refs=\"root\"/>\n";
	my @objects = @{$sets{$class}};
	print "<Objects objectIDs=\"" . join (" ", @objects) . "\"/>\n";
	print "</Classification>\n\n";
}


print "</ClassificationSet>

</ClassificationTree>\n";


exit(0);

