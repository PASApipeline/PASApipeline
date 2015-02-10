#!/usr/bin/env perl

use strict;

my $binparam;
unless ($binparam = $ARGV[0]) {
    die "usage: $0 bin_size [include zero counts]\n";
}

# for binning, want the following output:
#   range (0-10) -> 10
#   range (11-20) -> 20
#          (21-30) -> 30
#                         etc.
#  formula is:
#  x = number
#  binned(x) at seg-y = int ( (x-1+y)/y)

#analyse data

my $include_zero_counts = $ARGV[1];

my $value_counter = 0;
my %value_dist_counter;

my $entry_counter = 0;

my %tracker;
while (<STDIN>) {
    $_ =~ s/\s+//g;
    my $orig = $_;
    unless (/^[+-]?[\d\.]+$/) {
	print STDERR "can't process $_\n";
	next;
    }
    
    my $val_examine = $orig;
    if ($val_examine > 0) {
	$val_examine--;
    } elsif ($val_examine < 0) {
	$val_examine++;
    }
    
    my $index = int ( $val_examine/ $binparam); 
    if ($orig > 0) {
	$index++;
    } elsif ($orig < 0) {
	$index--;
    }
    
    $tracker{$index}++;

    $value_counter += abs($orig);
    $value_dist_counter{$index} += abs($orig);

    $entry_counter++;

}

my $total_counts = $entry_counter;

#print results
my @values = sort {$a<=>$b} keys %tracker;
my $maxvalue = $values[$#values];
my $minvalue = $values[0];

print STDERR "#bin\tcounts\t%counts\t%value of total\n";
for (my $i = $minvalue; $i <= $maxvalue; $i++) {
    my $index = ($i)  * ($binparam);
    my $value = (exists($tracker{$i})) ? $tracker{$i} : 0;
    
    if ($value != 0 || $include_zero_counts) {

	my $percentage_counts = sprintf ("%.2f", $value / $total_counts * 100);

	my $percentage_values = sprintf ("%.2f", $value_dist_counter{$i}/$value_counter * 100);
	print "$index\t$value\t$percentage_counts\t$percentage_values\n";
	
    }
}







