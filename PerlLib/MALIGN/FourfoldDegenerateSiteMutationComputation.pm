package MALIGN::FourfoldDegenerateSiteMutationComputation;

use strict;
use warnings;
use Carp;

use MALIGN::Aligned_sequence;

require Exporter;

our @ISA = qw (Exporter);
our @EXPORT = qw (compute4DsiteMutationRate);


my @FourD = qw( TC CT CC CG AC GT GC GG ); 
my %FDSites;
foreach my $dinuc (@FourD) {
	$FDSites{$dinuc}++;
}

####
sub compute4DsiteMutationRate {
	my ($alignedSeqA, $alignedSeqB) = @_;

	## note, these should be codon-aligned sequences!!!!

	my $seqA = uc $alignedSeqA->get_sequence();
	my $seqB = uc $alignedSeqB->get_sequence();
	
	if (length($seqA) != length($seqB)) {
		confess "error, seqsA and seqsB are of different lengths";
	}

	my $count_same = 0;
	my $count_diff = 0;

	my @codonsA = &_extract_codons($seqA);
	my @codonsB = &_extract_codons($seqB);

	for (my $i = 0; $i <= $#codonsA; $i++) {
		
		my $codonA = $codonsA[$i];
		my $codonB = $codonsB[$i];

		if ($codonA =~ /N/ || $codonB =~ /N/) {
			next;
		}

		my $codonA_prefix = substr($codonA, 0, 2);
		my $codonB_prefix = substr($codonB, 0, 2);

		if ($codonA_prefix eq $codonB_prefix && $FDSites{$codonA_prefix}) {
			
			if ($codonA eq $codonB) {
				$count_same++;
			}
			else {
				$count_diff++;
			}
		}
	}
	
	return ($count_same, $count_diff);
	
}


####
sub _extract_codons {
	my ($sequence) = @_;

	my $num_codons = length($sequence) / 3;

	my @codons;
	for (my $i = 0; $i < length($sequence); $i += 3) {
		my $codon = substr($sequence, $i, 3);
		push (@codons, $codon);
	}


	return (@codons);
}

1; #EOM

