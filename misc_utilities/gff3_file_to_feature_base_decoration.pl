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

my $usage = "\n\nusage: $0 gff3_file genome_db \n\n";

my $gff3_file = $ARGV[0] or die $usage;
my $fasta_db = $ARGV[1] or die $usage;

my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
    my $genome_seq = cdbyank_linear($asmbl_id, $fasta_db);

	my @pos_type_array = ();
	for (my $i = 0; $i < length($genome_seq); $i++) {
		$pos_type_array[$i] = 'G';
	}
    
    my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
	unless (@gene_ids) {
		die "Error, no gene_ids for [$asmbl_id] ";
	}

    foreach my $gene_id (@gene_ids) {
        my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
		
        my %params;
		
		#print $gene_obj_ref->toString();

	
				
        foreach my $isoform ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {

			#print "processing $isoform\n";
 
			my $orientation = $isoform->get_orientation();
			my ($model_lend, $model_rend) = sort {$a<=>$b} $isoform->get_model_span();
			my ($gene_lend, $gene_rend) = sort {$a<=>$b} $isoform->get_gene_span();
			
            my $isoform_id = $isoform->{Model_feat_name};
     

			my @exons = $isoform->get_exons();
			foreach my $exon (@exons) {
				my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
				for (my $i = $lend - 1; $i < $rend; $i++) {
					$pos_type_array[$i] = 'E';
				}
			}

			my @introns = $isoform->get_intron_coordinates();
			foreach my $intron (@introns) {
				my ($lend, $rend) = sort {$a<=>$b} @$intron;
				for (my $i = $lend; $i < $rend; $i++) {
					$pos_type_array[$i] = 'I';
				}
			}
			

        } # end of isoform

		

    } # end of gene

	print ">$asmbl_id\n" . join("", @pos_type_array) . "\n";


} # end of assembly


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


