#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 m2fmt_output  query_seq_lengths db_seq_lengths\n\n";

my $m2fmt_output = $ARGV[0] or die $usage;
my $query_seq_lengths = $ARGV[1] or die $usage;
my $db_seq_lengths = $ARGV[2] or die $usage;

main: {
	my %seq_lengths;
	{
		foreach my $seq_lengths_file ($query_seq_lengths, $db_seq_lengths) {
			open (my $fh, $seq_lengths_file) or die "Error, cannot open file $seq_lengths_file";
			while (<$fh>) {
				chomp;
				my ($seq_len, $acc) = split (/\t/);
				$acc =~ s/\s+$//;
				$seq_lengths{$acc} = $seq_len;
			}
			close $fh;
		}
	}
	
	open (my $fh, $m2fmt_output) or die "Error, cannot open file $m2fmt_output";
	while (<$fh>) {
		chomp;
		my @x = split (/\t/);

		## process query stats:
		my $query_acc = $x[0];
		my ($query_match_lend, $query_match_rend) = sort {$a<=>$b} ($x[17], $x[18]);
		my $query_match_len = $query_match_rend - $query_match_lend + 1;
		
		my $query_seq_len = $seq_lengths{$query_acc};
		my $percent_query_len = sprintf ("%.2f", $query_match_len / $query_seq_len * 100);

		push (@x, $query_seq_len, $percent_query_len);

		## process db stats:
		my $db_acc = $x[1];
		my ($db_match_lend, $db_match_rend) = sort {$a<=>$b} ($x[20], $x[21]);
		my $db_match_len = $db_match_rend - $db_match_lend + 1;
		
		my $db_seq_len = $seq_lengths{$db_acc};
		my $percent_db_len = sprintf ("%.2f", $db_match_len / $db_seq_len * 100);

		if ($query_acc eq $db_acc && $query_match_lend == $db_match_lend && $query_match_rend == $db_match_rend) {
			next; # self identical match, not relevant
		}
		
		push (@x, $db_seq_len, $percent_db_len);

		print join ("\t", @x) . "\n";
	}
	close $fh;
	
	
	exit(0);
	
}


		
			

