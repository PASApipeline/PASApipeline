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
use Overlap_piler;

my $usage = "\n\nusage: $0 gff3_file\n\n";

my $gff3_file = $ARGV[0] or die $usage;

my $gene_obj_indexer_href = {};

## Feature types wanted:  introns, exons, intergenes


my $introns_filename = "$gff3_file.introns";
open (my $introns_fh, ">$introns_filename") or die "Error, cannot write to $introns_filename";

my $exons_filename = "$gff3_file.exons";
open (my $exons_fh, ">$exons_filename") or die "Errror, cannot write to $exons_filename";

my $intergenes_filename = "$gff3_file.intergenes";
open (my $intergenes_fh, ">$intergenes_filename");

my $prime5utr_filename = "$gff3_file.5UTR";
open (my $prime5_fh, ">$prime5utr_filename") or die "Error, ccannot write to $prime5utr_filename";

my $prime3utr_filename = "$gff3_file.3UTR";
open (my $prime3_fh, ">$prime3utr_filename") or die "Error, cannot write to $prime3utr_filename";


## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
	my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
	my @gene_coords;


    foreach my $gene_id (@gene_ids) {
        my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
				
        foreach my $isoform ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
 
			my ($lend, $rend) = sort {$a<=>$b} $isoform->get_coords();
			push (@gene_coords, [$lend, $rend]);
			
			my $trans_id = $isoform->{Model_feat_name};
	
			print STDERR "// processing $asmbl_id, $gene_id, $trans_id\n";
			
		
			my $orient = $isoform->get_orientation();

			my @exons = $isoform->get_exons();
			
			foreach my $exon (@exons) {
				my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();

				print $exons_fh join("\t", $asmbl_id, $gene_id, $trans_id, $lend, $rend, $orient) . "\n";
			}
			
			
			my @introns = $isoform->get_intron_coordinates();
			
			foreach my $intron (@introns) {
				my ($lend, $rend) = sort {$a<=>$b} @$intron;
				
				print $introns_fh join("\t", $asmbl_id, $gene_id, $trans_id, $lend, $rend, $orient) . "\n";
			}


			my @prime5UTR = $isoform->get_5prime_UTR_coords();
			foreach my $utr (@prime5UTR) {
				my ($lend, $rend) = sort {$a<=>$b} @$utr;
				print $prime5_fh join("\t", $asmbl_id, $gene_id, $trans_id, $lend, $rend, $orient) . "\n";
			}
			
			my @prime3UTR = $isoform->get_3prime_UTR_coords();
			foreach my $utr (@prime3UTR) {
				my ($lend, $rend) = sort {$a<=>$b} @$utr;
				print $prime3_fh join("\t", $asmbl_id, $gene_id, $trans_id, $lend, $rend, $orient) . "\n";
			}
			
			
        }
    }
	
	
	my @gene_regions = &Overlap_piler::simple_coordsets_collapser(@gene_coords);
	
	@gene_regions = sort {$a->[0]<=>$b->[0]} @gene_regions;
	
	my $left_intergene = $gene_regions[0]->[0];
	foreach my $gene_region (@gene_regions) {
		my ($gene_lend, $gene_rend) = @$gene_region;
		
		if ($gene_lend > $left_intergene) {
			print $intergenes_fh join("\t", $asmbl_id, $left_intergene, $gene_lend - 1, "+") . "\n";
			print $intergenes_fh join("\t", $asmbl_id, $left_intergene, $gene_lend - 1, "-") . "\n";
		}
		
		$left_intergene = $gene_rend + 1;
	}
	
	
}

close $exons_fh;
close $introns_fh;
close $intergenes_fh;
close $prime5_fh;
close $prime3_fh;


exit(0);


