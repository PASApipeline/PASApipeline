#!/usr/bin/env perl

use strict;
use warnings;

use Gene_obj;

my $usage = "usage: $0 genes.gff\n\n";

my $gff = $ARGV[0] or die $usage;

my $counter = 0;

main: {

	open (my $fh, $gff) or die "Error, cannot open file $gff";
	while (<$fh>) {
		$counter++;
		
		chomp;
		my @x = split(/\t/);
	
		my $contig = $x[0];
		my ($lend, $rend) = ($x[3], $x[4]);
		my $orient = $x[6];
		my $info = $x[8];
		
		$info =~ /Name=(.*)/;
		my $name = $1 || "no name";
		
		my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);
		
		my $gene_obj = new Gene_obj();
		
		$gene_obj->{asmbl_id} = $contig;
		
		$gene_obj->{com_name} = $name;
		
		$gene_obj->{TU_feat_name} = "gene.$counter";
		$gene_obj->{Model_feat_name} = "model.$counter";
		
		my $coords_href = { $end5 => $end3 };
		
		$gene_obj->populate_gene_object($coords_href, $coords_href);
	
		print $gene_obj->to_GFF3_format() . "\n";
		
	}
	close $fh;



	exit(0);
}
