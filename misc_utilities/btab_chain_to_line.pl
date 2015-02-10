#!/usr/bin/env perl

use strict;
use warnings;

my @chain;
my $prev_chain;

while (<STDIN>) {
	chomp;
	my @x = split(/\t/);
	my ($contig, $acc, $con_end5, $con_end3, $acc_end5, $acc_end3, $per_id, $chain, $seg) = ($x[0], $x[5], $x[6], $x[7], $x[8], $x[9], $x[10], $x[13], $x[14]);
	
	if (defined($prev_chain) && $chain != $prev_chain) {
		
		&process_chain();
		@chain = ();
	}
	
	my $seg_struct = { contig => $contig,
					   acc => $acc,
					   con_end5 => $con_end5,
					   con_end3 => $con_end3,
					   acc_end5 => $acc_end5,
					   acc_end3 => $acc_end3,
					   chain => $chain,
					   seg => $seg,
					   perID => $per_id,
	};
	
	push (@chain, $seg_struct);
	
	$prev_chain = $chain;
}

&process_chain(); # get last one.

exit(0);

####
sub process_chain {
	
	my @contig_coords;
	my @seg_coords;
	my $contig;
	my $acc;
	my $sum_len_x_perid = 0;
	my $sum_len = 0;
	my $orient;
	
	foreach my $seg (@chain) {
		unless (defined $contig) {
			$contig = $seg->{contig};
			$acc = $seg->{acc};
		}
		my ($con_end5, $con_end3, $acc_end5, $acc_end3, $per_ID) = ($seg->{con_end5},
																	$seg->{con_end3},
																	$seg->{acc_end5},
																	$seg->{acc_end3},
																	$seg->{perID});
		
		unless ($orient) {
			if ($con_end5 != $con_end3 && $acc_end5 != $acc_end3) {
				my $genome_orient = ($con_end5 < $con_end3) ? "+" : "-";
				my $feature_orient = ($acc_end5 < $acc_end3) ? "+" : "-";
				
			    $orient = ($genome_orient eq $feature_orient) ? "+" : "-";
			}
		}

		my $len = abs($acc_end3 - $acc_end5) + 1;
		$sum_len_x_perid += $len*$per_ID;
		$sum_len += $len;
		push (@contig_coords, "$con_end5-$con_end3");
		push (@seg_coords, "$acc_end5-$acc_end3");
	}
	
	my $avg_perID =  $sum_len_x_perid / $sum_len;
	
	print "$contig\t$acc\t$orient\t" 
		. scalar(@contig_coords) . "\t" 
		. sprintf("%.3f", $avg_perID) . "\t"
		. join (",", @contig_coords) . "\t"
		. join (",", @seg_coords) 
		. "\n";
	
	return;
}
