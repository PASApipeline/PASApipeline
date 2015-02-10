#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 m2fmt_file max_stranded_tiers\n\n";

my $m2fmt_file = $ARGV[0] or die $usage;;
my $max_stranded_tiers = $ARGV[1] or die $usage;

unless ($max_stranded_tiers > 0) { die $usage; }

main: {
	my %contig_to_stranded_matches = &parse_matches($m2fmt_file);

	foreach my $contig_n_strand (keys %contig_to_stranded_matches) {
		my @matches = @{$contig_to_stranded_matches{$contig_n_strand}};

		@matches = reverse sort {$a->{bit_score}<=>$b->{bit_score}} @matches;

		## create the tiers list
		my @tiers;
		for (1..$max_stranded_tiers) {
			push (@tiers, []);
		}

		## add matches to tiers
		foreach my $match (@matches) {
		  tierSearch:
			foreach my $tier (@tiers) {
				if (! &has_overlapping_match($tier, $match)) {
					push (@$tier, $match);
					last tierSearch;
				}
			}
		}
		
		## report matches in tiers:
		foreach my $tier (@tiers) {
			foreach my $match (@$tier) {
				print $match->{line};
			}
		}
		
	}

	exit(0);
}

####
sub has_overlapping_match {
	my ($tier, $match) = @_;

	foreach my $other_match (@$tier) {
		my ($other_lend, $other_rend) = ($other_match->{lend}, $other_match->{rend});
		if ($match->{lend} <= $other_rend && $match->{rend} >= $other_lend) {
			return (1); # found overlap!
		}
	}

	return (0); # no overlapping match
}




####
sub parse_matches {
	my ($file) = @_;

	my %contig_to_matches;
	open (my $fh, $file) or die $!;
	while (<$fh>) {
		my $line = $_;
		my @x = split (/\t/);
		my ($contig, $bit_score, $lend, $rend, $orientA, $orientB) = ($x[0], $x[4], $x[17], $x[18], $x[16], $x[19]);

		($lend, $rend) = sort {$a<=>$b} ($lend, $rend);
		
		my $orient = ($orientA * $orientB > 0) ? '+' : '-';
		my $stranded_contig = "$contig$;$orient";

		push (@{$contig_to_matches{$stranded_contig}}, { bit_score => $bit_score,
														 lend => $lend,
														 rend => $rend,
														 line => $line,
													 } );
	}
	close $fh;

	return (%contig_to_matches);
}
