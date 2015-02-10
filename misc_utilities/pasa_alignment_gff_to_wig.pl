#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 pasa_alignment.gff3 [restrict_accs_file]\n\n";

my $pasa_gff = $ARGV[0] or die $usage;
my $restrict_accs_file = $ARGV[1];


my %restrict;

if ($restrict_accs_file) {
	my @accs = `cat $restrict_accs_file`;
	unless (scalar @accs > 1) {
		die $usage;
	}
	chomp @accs;
	
	%restrict = map {+ $_ => 1 } @accs;
}
	

my %chr_to_coverage;

open (my $fh, $pasa_gff) or die "Error, cannot open file $pasa_gff";
while (<$fh>) {
	chomp;
	my @x = split(/\t/);
	my $contig = $x[0];
	my $lend = $x[3];
	my $rend = $x[4];
	my $target_info = $x[8];
	
	$target_info =~ /Target=(asmbl_\d+)/;
	my $asmbl = $1 or die "Error, cannot extract pasa asmbl id from $target_info";
	
	if (%restrict) {
		unless ($restrict{$asmbl}) {
			next;
		}
	}

	
	
	for (my $i = $lend; $i <= $rend; $i++) {
		
		$chr_to_coverage{$contig}->[$i]++;
	}
	
	
}

print STDERR "\n";

## generate a tdf file
print "track name=\"Point coverage\" description=\"Sequence coverage\" visibility=2 itemRgb=\"On\"\n";

foreach my $chrom (sort keys %chr_to_coverage) {
	print "variableStep chrom=$chrom\n";
	
	my @pos_vals = @{$chr_to_coverage{$chrom}};
	
	for (my $i = 1; $i <= $#pos_vals; $i++) {
	   
		my $val = $pos_vals[$i] || 0;
		my $pos = $i;
		print "$pos\t$val\n";
		
	}
}

exit(0);


