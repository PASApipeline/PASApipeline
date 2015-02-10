#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;


my %genes;
my %gene_to_contig;

while (<STDIN>) {
	chomp;
	unless (/\w/) {next; }
	my @x = split (/\t/);
	my $contig = $x[0];
	my $source = $x[1];
	my $feat_type = $x[2];
	my $lend = $x[3];
	my $rend = $x[4];
	my $orient = $x[6];
	my $info = $x[8];
	
	

	my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);

	if ($feat_type =~ /initial|internal|terminal|single|stop_codon/) {
	
		$info =~ /gene_id \"([^\"]+)/ or die "Error, cannot parse gene_id from $info";
	
		my $gene_id = $1;
		

		$genes{$gene_id}->{$end5} = $end3;
		
		$gene_to_contig{$gene_id} = $contig;
	}



}

foreach my $gene_id (keys %genes) {
	my $coords_href = $genes{$gene_id};
	
	my $gene_obj = new Gene_obj();
	$gene_obj->populate_gene_object($coords_href, $coords_href);
	
	$gene_obj->{asmbl_id} = $gene_to_contig{$gene_id};
	$gene_obj->{TU_feat_name} = $gene_id;
	$gene_obj->{Model_feat_name} = "Model.$gene_id";
	
	$gene_obj->{com_name} = "predicted gene";
	
	$gene_obj->join_adjacent_exons();
	
	print $gene_obj->to_GFF3_format(source => "AUGUSTUS");
	print "\n";
}

exit(0);

