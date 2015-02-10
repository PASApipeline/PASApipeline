#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use CdbTools;
use GFF3_utils;
use Carp;
use Gene_validator;

my $usage = "\n\nusage: $0 gff3_file genome_db\n\n";

my $gff3_file = $ARGV[0] or die $usage;
my $fasta_db = $ARGV[1] or die $usage;


my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
    my $genome_seq = cdbyank_linear($asmbl_id, $fasta_db);
    
    my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
    foreach my $gene_id (@gene_ids) {
        my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
		
		my $errors = &validate_gene($gene_obj_ref, \$genome_seq);
		if ($errors) {
		  print "ERROR: $gene_id\t$errors\n";
		}
		else {
		  print "OK: $gene_id\n";
		}
		
    }
}


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


