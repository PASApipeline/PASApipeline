#!/usr/bin/env perl

use strict;
use warnings;

use PSL_parser;
use Fasta_reader;

use Overlap_piler;

my $usage = "usage: $0 file.psl.top target_seqs.fa\n\n";

my $psl_file = $ARGV[0] or die $usage;
my $target_seqs_fa = $ARGV[1] or die $usage;

main: {

	
	my %acc_to_spans;

	my %seq_lengths = &parse_seq_lengths($target_seqs_fa);

	my $psl_parser = new PSL_parser($psl_file);

	while (my $psl_entry = $psl_parser->get_next()) {

		my $per_id = $psl_entry->get_per_id();
		
		unless ($per_id >= 95) { next; }
		
		my $cds_name = $psl_entry->get_T_name();
		
		my ($g_coords_aref, $q_coords_aref) = $psl_entry->get_alignment_coords();

		foreach my $coordset (@$g_coords_aref) {
			
			push (@{$acc_to_spans{$cds_name}}, $coordset);
		}

	}

	foreach my $acc (keys %seq_lengths) {

		my $target_length = $seq_lengths{$acc};
		
		
		if (my $aref = $acc_to_spans{$acc}) {

			my @piles = &Overlap_piler::simple_coordsets_collapser(@$aref);

			my $sum_len = 0;

			foreach my $pile (@piles) {
				my ($lend, $rend) = @$pile;
				
				$sum_len += $rend - $lend + 1;
				
			}

			my $percent_coverage = sprintf("%.2f", $sum_len / $target_length * 100);

			print "$acc\t$percent_coverage\n";
		}
		else {
			print "$acc\t0\n";
		}
	}
	

	exit(0);
	
	
	

}


####
sub parse_seq_lengths {
	my ($fasta) = @_;

	
	my %lengths;
	

	my $fasta_reader = new Fasta_reader($fasta);

	while (my $seq_obj = $fasta_reader->next()) {

		my $acc = $seq_obj->get_accession();

		my $seq = $seq_obj->get_sequence();

		$lengths{$acc} = length($seq);
	}

	return(%lengths);
}


