#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");

use Gene_obj;
use Fasta_reader;

use Data::Dumper;

my $usage = "usage: $0 annots.bed genome.fasta\n\n";

my $annots_bed = $ARGV[0] or die $usage;
my $genome_fasta = $ARGV[1] or die $usage;


my $DEBUG = 0;

main: {
	
	my $fasta_reader = new Fasta_reader($genome_fasta);

	my %genome = $fasta_reader->retrieve_all_seqs_hash();

	my $counter = 0;

	open (my $fh, $annots_bed) or die "Error, cannot open file $annots_bed";
	while (<$fh>) {
		my $line = $_;
		print if $DEBUG;
		
		my @x = split(/\t/);

		my $scaff = $x[0];
		my $gene_lend = $x[1] + 1;
		my $gene_rend = $x[2];

		my $com_name = $x[3];

		my $score = $x[4];
		my $orient = $x[5];
		
		if ($orient eq '*') {
			$orient = '+';
		}
		

		my $coding_lend = $x[6] + 1;
		my $coding_rend = $x[7];

		my $rgb_color = $x[8];

		my $num_exons = $x[9];

		my $lengths_text = $x[10];
		my $exon_relative_starts_text = $x[11];

		my @lengths = split(/,/, $lengths_text);
		my @exon_relative_starts = split(/,/, $exon_relative_starts_text);

		my @exons;

		while (@lengths) {
			my $len = shift @lengths;
			my $start = shift @exon_relative_starts;
			
			my $exon_lend = $gene_lend + $start;
			my $exon_rend = $exon_lend + $len - 1;
			

			print "Len: $len, start=$start   ====>  $exon_lend - $exon_rend\n" if $DEBUG;

			push (@exons, [$exon_lend, $exon_rend]);
			
		}

		
		print "Coding: $coding_lend-$coding_rend, Exons: " . Dumper (\@exons) if $DEBUG;
		

		eval {
			my $gene_obj = new Gene_obj();
			$gene_obj->build_gene_obj_exons_n_cds_range(\@exons, $coding_lend, $coding_rend, $orient);
			
			$gene_obj->{com_name} = $com_name;
			$gene_obj->{asmbl_id} = $scaff;
			
			$counter++;
			$gene_obj->{TU_feat_name} = "gene.$counter";
			$gene_obj->{Model_feat_name} = "model.$counter";
			
			print $gene_obj->to_GFF3_format() . "\n" if $DEBUG;
			
			my $cdna_seq = $gene_obj->create_cDNA_sequence(\$genome{$scaff});
			
			$com_name =~ s/\s+/_/g;
			
			print ">$com_name\n$cdna_seq\n";
			
		};

		if ($@) {
			print STDERR "Error, couldn't reconstruct sequence for entry: $line\n";
		}
		
	}
	
	exit(0);
}
	
