#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use PSL_parser;
use Gene_obj;

my $usage = "usage: $0 file.psl(x)\n\n";

my $psl_file = $ARGV[0] or die $usage;

main: {

	my $psl_parser = new PSL_parser($psl_file);
	
	while (my $pe = $psl_parser->get_next()) {

		my $genome_contig = $pe->get_T_name();
		my $query_name = $pe->get_Q_name();

		my $strand = $pe->get_strand();
		
		my ($genome_coords_aref, $query_coords_aref) = $pe->get_alignment_coords();

		my %coords;
		
		foreach my $coords_aref (@$genome_coords_aref) {
			
			my ($lend, $rend) = @$coords_aref;

			my ($end5, $end3) = ($strand eq '+') ? ($lend, $rend) : ($rend, $lend);

			$coords{$end5} = $end3;
		}

		my $gene_obj = new Gene_obj();
		$gene_obj->populate_gene_object(\%coords, \%coords);
		$gene_obj->{asmbl_id} = $genome_contig;
		$gene_obj->{com_name} = $query_name;

		print $gene_obj->to_BED_format();

	}

	exit(0);
}


		
