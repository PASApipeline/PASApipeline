#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use CdbTools;
use Nuc_translator;

my $usage = "\nusage: $0 range_file genomeMultiFasta [flank=0]\n\n"
	. "format for range_file:\n"
	. "accession (tab) end5 (tab) end3\n\n\n";

my $range_file = $ARGV[0] or die $usage;
my $genome = $ARGV[1] or die $usage;
my $flank = $ARGV[2] || 0;


main: {
	my %genome_acc_to_coord_sets;

	open (my $fh, $range_file) or die "Error, cannot open file $range_file";
	while (<$fh>) {
		chomp;
		my ($acc, $end5, $end3) = split (/\t/);

		push (@{$genome_acc_to_coord_sets{$acc}}, [$end5, $end3]);
	}
	close $fh;

	foreach my $acc (sort keys %genome_acc_to_coord_sets) {
		my $feature_list_aref = $genome_acc_to_coord_sets{$acc};

		my $genome_seq = &cdbyank_linear($acc, $genome);

		foreach my $coordset (sort {$a->[0]<=>$b->[0]} @$feature_list_aref) {
			my ($end5, $end3) = @$coordset;
			
			my ($lend, $rend) = sort {$a<=>$b} ($end5, $end3);
			
			$lend -= $flank;
			$rend += $flank;

			$lend = 1 if $lend < 1;
			

			my $orient = ($end5 < $end3) ? '+' : '-';
			
			my $subseq = substr($genome_seq, $lend - 1, $rend - $lend + 1);
			
			if ($orient eq '-') {
				$subseq = &reverse_complement($subseq);
			}

			$subseq =~ s/(\S{60})/$1\n/g;
			chomp $subseq;
			
			print ">$acc-$end5-$end3 flank:$flank\n$subseq\n";
					
		}
	}

	exit(0);
}
