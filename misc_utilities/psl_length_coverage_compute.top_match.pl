#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});

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

	my %gene_to_hit_list;

	while (my $psl_entry = $psl_parser->get_next()) {

		my $per_id = $psl_entry->get_per_id();
		
		unless ($per_id >= 95) { next; }
		
		my $cds_name = $psl_entry->get_T_name();

		my $score = $psl_entry->get_match_count() * 5 - $psl_entry->get_mismatch_count() * 4;
		
		$psl_entry->{_score} = $score;

		push (@{$gene_to_hit_list{$cds_name}}, $psl_entry);
		
	}

	
	foreach my $cds (keys %seq_lengths) {
		
		if (exists $gene_to_hit_list{$cds}) {
			
			my @psl_entries = @{$gene_to_hit_list{$cds}};
			
			@psl_entries = sort {$a->{_score}<=>$b->{_score}} @psl_entries;
			
			my $target_length = $seq_lengths{$cds};
			
			my $sum_len = 0;
			
			my $top_psl = pop @psl_entries;
			
			my ($genome_aref, $trans_aref) = $top_psl->get_alignment_coords();
			
			
			foreach my $segment (@$genome_aref) {
				my ($lend, $rend) = @$segment;
				
				$sum_len += ($rend - $lend) + 1;
				
			}
			
			my $percent_coverage = sprintf("%.2f", $sum_len / $target_length * 100);
			
			print "$cds\t$percent_coverage\n";
		}
		
		else {
			print "$cds\t0\n";
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


