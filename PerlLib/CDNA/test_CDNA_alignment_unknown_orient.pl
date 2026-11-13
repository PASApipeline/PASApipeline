#!/usr/bin/env perl

# Regression test for issue #345: an alignment with unknown spliced orientation ('?')
# must be built in both orientations, keeping the one with the longer ORF.
# Adapted from the reproducer contributed by JMusser499 in #345.
#
# usage: perl PerlLib/CDNA/test_CDNA_alignment_unknown_orient.pl    (exit 0 = pass)

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/..");
use CDNA::CDNA_alignment;
use CDNA::Alignment_segment;

{   # records the orientations tried, then runs PASA's own implementation unchanged
    package SyntheticAlignment;
    our @ISA = ('CDNA::CDNA_alignment');
    our @calls;
    sub _get_gene_obj_via_alignment_by_orient {
        my ($self, $orient, $partial) = @_;
        push @calls, $orient;
        return $self->SUPER::_get_gene_obj_via_alignment_by_orient($orient, $partial);
    }
}

# 60 nt, one exon aligned on the minus strand:
#   forward ORF ATG-(GCC)x3-TAA (15 nt), reverse ORF ATG-(AAA)x9-TAA (33 nt)
my $dna = 'CCC' . 'ATG' . ('GCC' x 3) . 'TAA' . 'CCCCCC'
        . 'TTA' . ('TTT' x 9) . 'CAT' . 'CCC';

my %expected = (
    '+' => ['+', '1-60', '4-18', 15, 'MAAA*', '+'],
    '-' => ['-', '1-60', '57-25', 33, 'MKKKKKKKKK*', '-'],
    '?' => ['-', '1-60', '57-25', 33, 'MKKKKKKKKK*', '+,-'],
);
my @fields = ('gene strand', 'exon interval', 'CDS 5prime-3prime', 'CDS nt', 'protein', 'orientations tried');

my $failures = 0;
foreach my $orient ('+', '-', '?') {
    my $segment = CDNA::Alignment_segment->new(length($dna), 1, 1, length($dna), 100);
    my $alignment = SyntheticAlignment->new(length($dna), [$segment], \$dna);
    die "fixture is not a minus-aligned single exon of unknown spliced orientation\n"
        unless $alignment->get_aligned_orientation() eq '-' && $alignment->get_spliced_orientation() eq '?';
    $alignment->set_spliced_orientation($orient);
    @SyntheticAlignment::calls = ();

    my $gene = $alignment->get_gene_obj_via_alignment();
    my @observed = ($gene->get_orientation(),
                    join(',', map { join('-', sort {$a <=> $b} $_->get_coords()) } $gene->get_exons()),
                    join('-', $gene->get_CDS_span()),
                    $gene->get_CDS_length(),
                    $gene->get_protein_sequence(),
                    join(',', @SyntheticAlignment::calls));

    foreach my $i (0..$#fields) {
        my $ok = $observed[$i] eq $expected{$orient}->[$i];
        print join("\t", $orient, $fields[$i], $expected{$orient}->[$i], $observed[$i], $ok ? "ok" : "FAIL") . "\n";
        $failures++ unless $ok;
    }
}

print "failures: $failures\n";
exit($failures ? 1 : 0);
