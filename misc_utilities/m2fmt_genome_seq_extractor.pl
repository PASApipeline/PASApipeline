#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use CdbTools;
use Nuc_translator;

my $usage = "\nusage: $0 genomeDB m2fmt_file [leftflank:rightflank]\n\n";

my $genomeDB = $ARGV[0] or die $usage;
my $m2fmt = $ARGV[1] or die $usage;
my $flank = $ARGV[2];

my ($left_flank, $right_flank) = (0,0);

if ($flank) {
    if ($flank =~ /^(\d+)\:(\d+)$/) {
		($left_flank, $right_flank) = ($1, $2);
	}
	else {
		die $usage;
	}
}

my %contig_to_structs;

open (my $fh, $m2fmt) or die "Error, cannot read $m2fmt";
while (<$fh>) {
	chomp;
	my @x = split (/\t/);
	my ($acc, $contig, $acc_lend, $acc_rend, $contig_end5, $contig_end3) = ($x[0], $x[1], $x[17], $x[18], $x[20], $x[21]);
	
	my ($acc_orient, $contig_orient) = ($x[16], $x[19]);
	
	push (@{$contig_to_structs{$contig}}, { acc => $acc,
											acc_lend => $acc_lend,
											acc_rend => $acc_rend,
											acc_orient => $acc_orient,
											
											contig_end5 => $contig_end5,
											contig_end3 => $contig_end3, 
											contig_orient => $contig_orient,
										} );
	
}
close $fh;


my $counter = 0;

foreach my $contig (keys %contig_to_structs) {
	
	my $genome_seq = &cdbyank_linear($contig, $genomeDB);
	
	foreach my $struct (@{$contig_to_structs{$contig}}) {

		my ($acc, $acc_lend, $acc_rend, $contig_end5, $contig_end3) = ($struct->{acc},
																	   $struct->{acc_lend},
																	   $struct->{acc_rend},
																	   $struct->{contig_end5},
																	   $struct->{contig_end3});

		my ($acc_orient, $contig_orient) = ($struct->{acc_orient}, $struct->{contig_orient});
		

		my $orient = ($acc_orient * $contig_orient > 0) ? '+' : '-';
		
		my ($contig_lend, $contig_rend) = sort {$a<=>$b} ($contig_end5, $contig_end3);

		if ($flank) {
			if ($orient eq '+') {
				$contig_lend -= $left_flank;
				$contig_rend += $right_flank;
			}
			else {
				$contig_lend -= $right_flank;
				$contig_rend += $left_flank;
			}
			
			if ($contig_lend < 1) {
				$contig_lend = 1;
			}

		}
		

		my $geneseq = substr($genome_seq, $contig_lend -1, $contig_rend - $contig_lend + 1);

		if ($orient eq '-') {
			$geneseq = &reverse_complement($geneseq);
		}
		
		$counter++;

		print ">s$counter $acc $acc_lend-$acc_rend $contig:$contig_end5-$contig_end3\n$geneseq\n";
	}
}


exit(0);


