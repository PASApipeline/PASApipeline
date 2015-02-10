#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;

use List::Util qw (shuffle);

my $usage = "usage: $0 CDS hexamerScores [randIter=100]\n\n";

my $cds_file = $ARGV[0] or die $usage;
my $hexamer_scores_file = $ARGV[1] or die $usage;
my $NUM_ITER = $ARGV[2] || 100;

main: {
	
	my %scores = &parse_hexamer_scores($hexamer_scores_file);
	
	my $fasta_reader = new Fasta_reader($cds_file);
	while (my $seq_obj = $fasta_reader->next()) {
		
		my $accession = $seq_obj->get_accession();
		my $sequence = uc $seq_obj->get_sequence();

		my $seq_length = length($sequence);
		
		if ($seq_length < 5) {
			next;
		}
		
		my $score = &score_sequence($sequence, \%scores);

		my @rand_scores;
		for my $it (1..$NUM_ITER) {
			#print STDERR "\riter[$it]   ";
			$sequence = join("", shuffle( split(//, $sequence) ) );

			push (@rand_scores, &score_sequence($sequence, \%scores));
		}
		
		my @rand_scores_above = grep { $_ >= $score } @rand_scores;

		my $p_value = scalar(@rand_scores_above) / $NUM_ITER;

		print "$accession\t" . sprintf("%.2f", $score) . "\t" . sprintf("%.3f", $p_value) . "\n";
		#print "\t" . join(", ", @rand_scores) . "\n";
	}
	
}

####
sub score_sequence {
	my ($sequence, $scores_href) = @_;
	
	my $seq_length = length($sequence);
	
	## init score to first pentamer
	my $pentamer = substr($sequence, 0, 5);
	my $framed_pentamer = "${pentamer}-0";
	my $score = $scores_href->{$framed_pentamer} || 0;
	
	for (my $i = 5; $i <= $seq_length - 6; $i++) {
		my $hexamer = substr($sequence, $i, 6);
		my $frame = $i % 3;
		my $framed_hexamer = "${hexamer}-${frame}";
		my $hex_score = $scores_href->{$framed_hexamer} || 0;
		$score += $hex_score;
	}
	
	return($score);
}


####
sub parse_hexamer_scores {
	my ($hexamer_scores_file) = @_;

	my %scores;
	open (my $fh, $hexamer_scores_file) or die "Error, cannot open $hexamer_scores_file";
	while (<$fh>) {
		chomp;
		my ($token, $score) = split (/\t/);
		$scores{$token} = $score;
	}
	close $fh;

	return (%scores);
}

