#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Gene_obj;

my $usage = "usage: $0 m2fmt_file\n\n";

my $m2fmt_file = $ARGV[0] or die $usage;

my %seen;
open (my $fh, $m2fmt_file) or die "Error, cannot open file $m2fmt_file";
while (<$fh>) {
	chomp;
	my @x = split (/\t/);
	my ($acc, $contig_acc, $match_lend, $match_rend, $contig_lend, $contig_rend) = ($x[0], $x[1], $x[17], $x[18], $x[20], $x[21]);
	
	my ($end5, $end3) = ($match_lend < $match_rend) ? ($contig_lend, $contig_rend) : ($contig_rend, $contig_lend);
	
	if ($seen{$acc}) {
		$acc .= "." . ++$seen{$acc};
	}
	else {
		$seen{$acc}++;
	}
	
	my $gene_obj = new Gene_obj();
	$gene_obj->populate_gene_object({$end5=>$end3}, {$end5=>$end3});

	$gene_obj->{asmbl_id} = $contig_acc;
	$gene_obj->{TU_feat_name} = $acc;
	$gene_obj->{Model_feat_name} = "$acc.mRNA";
	$gene_obj->{com_name} = "";
	
	print $gene_obj->to_GFF3_format() . "\n";

}

exit(0);


