#!/usr/bin/env perl

use strict;
use warnings;

# calculates a Jaccard similarity coefficient based on unique 7-mer composition of query and database seqs.

use lib ($ENV{EUK_MODULES});
use Fasta_reader;

my $KMER_SIZE = 7;


my $usage = "usage: $0 database query [numTopHits=10]\n\n";

my $databaseFile = $ARGV[0] or die $usage;
my $queryFile = $ARGV[1] or die $usage;
my $num_top_hits = $ARGV[2] || 10;

main: {
		
	my @scores;

	my $query_reader = new Fasta_reader($queryFile);
	my $querySeqObj = $query_reader->next();
	my $query_seq = uc $querySeqObj->get_sequence();
	my $query_acc = $querySeqObj->get_accession();
	
	$query_seq =~ s/U/T/g; 
	
	my %unique_query_Kmers = &get_unique_Kmers($query_seq);
	
	my $counter = 0;
	
	my $db_reader = new Fasta_reader($databaseFile);
	while (my $seqObj = $db_reader->next()) {
	
		$counter++;
		print STDERR "\ranalyzing seq $counter          ";
		
		my $acc = $seqObj->get_accession();
		my $sequence = uc $seqObj->get_sequence();
		$sequence =~ s/U/T/g;
		
		my %unique_db_Kmers = &get_unique_Kmers($sequence);
		
		my $similarity_coeff = &compute_similarity_coeff(\%unique_query_Kmers, \%unique_db_Kmers);
		
		if (scalar(@scores) < $num_top_hits) {
			push (@scores, [$acc, $similarity_coeff]);
			@scores = reverse sort {$a->[1]<=>$b->[1]} @scores;
		}
		elsif ($scores[$#scores]->[1] < $similarity_coeff) {
			pop @scores; # remove lowest scoring entry
			push (@scores, [$acc, $similarity_coeff]);
			@scores = reverse sort {$a->[1]<=>$b->[1]} @scores;
		}
	}
	
	## refresh terminal:
	print STDERR "\r                                                                        \r";
	
	@scores = reverse sort {$a->[1]<=>$b->[1]} @scores;
	
	for (my $i = 0; $i < $num_top_hits && $i <= $#scores; $i++) {
		my $score_struct = $scores[$i];
		my ($acc, $score_val) = @$score_struct;
		
		$score_val = sprintf("%.5f", $score_val);
		print "$score_val\t$acc\t$query_acc\n";
	}
	
	exit(0);
}



####
sub compute_similarity_coeff {
	my ($kmers_A_href, $kmers_B_href) = @_;
	
	my %all_kmers;
	foreach my $kmer (keys %$kmers_A_href, keys %$kmers_B_href) {
		$all_kmers{$kmer}++;
	}
	
	my $kmer_intersection = 0;
	my @all_kmers_list = keys %all_kmers;
	my $num_all_kmers = scalar (@all_kmers_list);
	
	foreach my $kmer (@all_kmers_list) {
		if ($all_kmers{$kmer} > 1) {
			$kmer_intersection++;
		}
	}
	
	my $similarity_coeff = $kmer_intersection / $num_all_kmers;
	
	return ($similarity_coeff);
}


####
sub get_unique_Kmers {
	my ($sequence) = @_;

	# print "SEQ: $sequence\n";
	
	my %kmer_counts;
	my @chars = split (//, $sequence);
	for (my $i = 0; $i + $KMER_SIZE <= length($sequence); $i++) {
		my $word = join ("", @chars[$i..$i+$KMER_SIZE-1]);
		$kmer_counts{$word}++;
		# print "$word\n";
	}
	
	## remove those that occur multiple times.
	foreach my $kmer (keys %kmer_counts) {
		if ($kmer_counts{$kmer} > 1) {
			delete($kmer_counts{$kmer});
		}
	}

	return (%kmer_counts);
}
	
	
	
	
