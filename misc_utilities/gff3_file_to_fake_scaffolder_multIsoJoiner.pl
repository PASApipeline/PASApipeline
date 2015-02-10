#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use CdbTools;
use GFF3_utils;
use Carp;
use Nuc_translator;

my $usage = "\n\nusage: $0 gff3_file genome_db [N_sep=100]\n\n";

my $gff3_file = $ARGV[0] or die $usage;
my $fasta_db = $ARGV[1] or die $usage;

my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

my $fake_contig = "";


foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
    my $genome_seq = cdbyank_linear($asmbl_id, $fasta_db);
    
    my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
    foreach my $gene_id (@gene_ids) {
        my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
		
        my %exon_coords;
		
		my $isoform_names = "";
	
        foreach my $isoform ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
			
			if ($isoform_names) {
				$isoform_names .= ";";
			}
			
			$isoform_names .= $isoform->{Model_feat_name} . " ";
			
			my @exons = $isoform->get_exons();
			
			foreach my $exon (@exons) {
			
				my ($end5, $end3) = $exon->get_coords();
				$exon_coords{$end5}=$end3;
			}

		}

		my $new_gene_obj = new Gene_obj();
		$new_gene_obj->populate_gene_object(\%exon_coords, \%exon_coords);
		
		$new_gene_obj->join_adjacent_exons();

		my $fake_cds = $new_gene_obj->create_CDS_sequence(\$genome_seq);
				
		if ($gene_obj_ref->get_orientation() eq '-') {
			$fake_cds = &reverse_complement($fake_cds);
		}

		if ($fake_contig) {
			$fake_contig .= 'N' x 100;
		}
		my $start = length($fake_contig) + 1;
		$fake_contig .= $fake_cds;

		my $end = length($fake_contig);
		
		
		
		my $name = $gene_id . "; $isoform_names;" . $gene_obj_ref->{com_name};
		print join("\t", "fake", ".", "gene", $start, $end, ".", '+', '.', 
				   $name) . "\n";
	}
	
}

open (my $ofh, ">fake.$$.contig") or die $!;
$fake_contig =~ s/(\S{60})/$1\n/g;
print $ofh ">fake\n$fake_contig\n";
close $ofh;


exit(0);


####
sub add_flank {
	my ($seq, $upstream_flank, $downstream_flank, $lend, $rend, $orientation, $genome_seq_ref) = @_;
	
	my $far_left = ($orientation eq '+') ? $lend - $upstream_flank : $lend - $downstream_flank;
	
	if ($far_left < 1) { $far_left = 1; }
	
	my $flank_right = ($orientation eq '+') ? $downstream_flank : $upstream_flank;

	my $left_seq = substr($$genome_seq_ref, $far_left - 1, $lend - $far_left);

	my $right_seq = substr($$genome_seq_ref, $rend, $flank_right);
	
	if ($orientation eq '+') {
		return (lc($left_seq) . uc($seq) . lc($right_seq));
	}
	else {
		return (lc(&reverse_complement($right_seq)) . uc($seq) . lc(&reverse_complement($left_seq)));
	}
}


