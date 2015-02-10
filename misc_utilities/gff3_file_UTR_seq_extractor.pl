#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use CdbTools;
use GFF3_utils;
use Carp;
use Fasta_reader;

my $usage = "\n\nusage: $0 gff3_file genome.fasta primeUTR[5|3]\n\n";

my $gff3_file = $ARGV[0] or die $usage;
my $genome_fasta = $ARGV[1] or die $usage;
my $primeUTR = $ARGV[2] or die $usage;

unless ($primeUTR == 5 || $primeUTR==3) {
	die $usage;
}

my $fasta_reader = new Fasta_reader($genome_fasta);

my %genome_seqs = $fasta_reader->retrieve_all_seqs_hash();

my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
    my $genome_seq = $genome_seqs{$asmbl_id};

    my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
    foreach my $gene_id (@gene_ids) {
        my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};

		my $counter = 0;
		
		foreach my $isoform ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
			
			my $TU = $isoform->{TU_feat_name};
			my $model = $isoform->{Model_feat_name};
			
			my $utr_seq = ($primeUTR == 5) ? $isoform->get_5prime_UTR_sequence(\$genome_seq)
				: $isoform->get_3prime_UTR_sequence(\$genome_seq);
			
			if ($utr_seq) {
				print ">$model-${primeUTR}utr\t$model\t$TU\t${primeUTR}-primeUTR\n$utr_seq\n";
			}
		}
			
	}
		
				
	
}



exit(0);

