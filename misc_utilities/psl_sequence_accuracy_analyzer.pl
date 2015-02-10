#!/usr/bin/env perl

use strict;
use warnings;

use PSL_parser;

my $usage = "usage: $0 blat.pslx.top MIN_INTRON_LENGTH\n\n";

my $psl_file = $ARGV[0] or die $usage;
my $min_intron_length = $ARGV[1] or die $usage;

main: {

	
	my $psl_reader = new PSL_parser($psl_file);

	my $sum_alignments = 0;
	my $sum_matches = 0;
	my $sum_mismatches = 0;

	my $sum_genome_insertions = 0;
	my $sum_transcript_insertions = 0;

	while (my $entry = $psl_reader->get_next()) {

		$sum_alignments++;

		my $mismatch_count = $entry->get_mismatch_count();

		$sum_mismatches += $mismatch_count;
		
		my $match_count = $entry->get_match_count();
		
		$sum_matches += $match_count;
		
		
		my ($genome_coords_aref, $trans_coords_aref) = $entry->get_alignment_coords();

		my @genome_deltas = &get_deltas($genome_coords_aref);
		
		my @trans_deltas = &get_deltas($trans_coords_aref);
	
		foreach my $genome_delta (@genome_deltas) {
			if ($genome_delta < $min_intron_length) {
				$sum_genome_insertions += $genome_delta;
			}
		}

		foreach my $trans_delta (@trans_deltas) {
			$sum_transcript_insertions += $trans_delta;
		}

	}
	

	my $total_aligned_bases = $sum_matches + $sum_mismatches;

	print "#transcripts\taligned_bases\tmatches\tmismatches\tgenome_insert\ttranscript_insert\n";
	print "$sum_alignments\t$total_aligned_bases\t$sum_matches\t$sum_mismatches\t$sum_genome_insertions\t$sum_transcript_insertions\n\n";	
	
	my $mismatch_rate = sprintf("%.2e", $sum_mismatches / $total_aligned_bases);
	my $genome_insert_rate = sprintf("%.2e", $sum_genome_insertions / $total_aligned_bases);
	my $transcript_insert_rate = sprintf("%.2e", $sum_transcript_insertions / $total_aligned_bases);
	
	print "\nMismatch rate: $mismatch_rate\n"
		. "Genome insertion rate: $genome_insert_rate\n"
		. "Transcript insert rate: $transcript_insert_rate\n\n";
	
	exit(0);

}


####
sub get_deltas {
	my ($coords_aref) = @_;

	my @coordsets = @$coords_aref;
	if (scalar @coordsets == 1) {
		return();
	}

	foreach my $coordset (@coordsets) {
		@$coordset = sort {$a<=>$b} @$coordset;
	}

	@coordsets = sort {$a->[0]<=>$b->[0]} @coordsets;

	my @deltas;

	for (my $i = 1; $i <= $#coordsets; $i++) {

		my $prev_coordset = $coordsets[$i-1];
		
		my $curr_coordset = $coordsets[$i];

		my $delta = $curr_coordset->[0] - $prev_coordset->[1] - 1;

		push (@deltas, $delta) if $delta;
	}

	return(@deltas);
}
			
