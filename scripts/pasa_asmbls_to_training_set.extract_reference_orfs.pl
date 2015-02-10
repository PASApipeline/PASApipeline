#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use GFF3_utils;
use TierFeatures;

my $usage = "usage: $0 best_candidates.gff3 [minProtLength=300]\n\n";

my $gff3_file = $ARGV[0] or die $usage;
my $min_prot_length = $ARGV[1] || 300;

my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
    my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
	my @features;

    foreach my $gene_id (@gene_ids) {
        my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
                
		my ($lend, $rend) = sort {$a<=>$b} $gene_obj_ref->get_coords();
		my $cds_length = $gene_obj_ref->get_CDS_length();
		my $prot_length = $cds_length / 3;

		unless ($prot_length > $min_prot_length) { next; }

		my $feature = new TierFeatures::Feature($lend, $rend, $gene_id);

		$feature->{__PROT_LENGTH} = $prot_length;

		push (@features, $feature);
	}

	@features = reverse sort {$a->{__PROT_LENGTH}<=>$b->{__PROT_LENGTH}} @features;

	my $FT = new TierFeatures();
	my @tiers = $FT->tier_features(@features);

	my $first_tier = $tiers[0];
	
	foreach my $ele (@$first_tier) {

		my $gene_id = $ele->{feat_id};

		my $gene_obj = $gene_obj_indexer_href->{$gene_id};
		
		print $gene_obj->to_GFF3_format() . "\n";
	}

}


exit(0);


		
