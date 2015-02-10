#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 blast.m8\n\n";

my $m8 = $ARGV[0] or die $usage;


my %best_hits;

open (my $fh, $m8) or die "Error, cannot open file $m8";
while (<$fh>) {
	chomp;
	my @x = split(/\t/);

	my $accA = $x[0];
	my $accB = $x[1];
	my $per_id = $x[2];
	my $e_value = $x[10];
	my $blast_score = $x[11];
	
	if ($accA eq $accB) { next; }
	
	my $struct = $best_hits{$accA};
		
	if ( (! $struct) || $struct->{score} < $blast_score) {
		
		$best_hits{$accA} = { accA => $accA,
							  accB => $accB,
							  score => $blast_score,
							  per_id => $per_id,
							  e_value => $e_value,
						  };
		
		

	}
}

close $fh;

my %seen;

use Data::Dumper;

foreach my $struct (values %best_hits) {
	
	my $acc_A = $struct->{accA};
	my $acc_B = $struct->{accB};
	my $score = $struct->{score};
	my $per_id = $struct->{per_id};
	my $e_value = $struct->{e_value};
	
	if ($seen{$acc_A} || $seen{$acc_B}) { next; }
	
	if (my $other_struct = $best_hits{$acc_B}) {
		
		my $best_match = $other_struct->{accB};

		#print "$acc_A\t$acc_B\t$best_match\n";
		
		if ($other_struct->{accB} eq $acc_A) {

			print join("\t", $acc_A, $acc_B, $per_id, $e_value, $score) . "\n";
			
			$seen{$acc_A} = 1;
			$seen{$acc_B} = 1;

		}
	}


}


exit(0);


		
