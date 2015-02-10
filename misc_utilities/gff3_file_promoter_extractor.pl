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
use Nuc_translator;
use Getopt::Long qw(:config no_ignore_case bundling);


my $usage = <<_EOUSAGE_;

##############################################################
#
#  Required:
#
#  --genes          gene annotations in GFF3 format
#  --genome         genome sequence in multi-FASTA format
#
# Optional:
#
#  --require_5prime_UTR   only extracts sequences corresponding to genes with 5' UTRs
#
#  --upstream            default: 500
#  --downstream          default: 100
#
###############################################################

_EOUSAGE_

	;


my $gff3_file;
my $genome_fasta;
my $require_5prime_UTR_flag = 0;
my $upstream = 500;
my $downstream = 100;


&GetOptions( 'genes=s' => \$gff3_file,
			 'genome=s' => \$genome_fasta,
			 'require_5prime_UTR_flag' => \$require_5prime_UTR_flag,
			 'upstream=i' => \$upstream,
			 'downstream=i' => \$downstream,
			 
			 );


unless ($gff3_file && $genome_fasta) {
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
		

		my $TU = $gene_obj_ref->{TU_feat_name};
		
		my $model = $gene_obj_ref->{Model_feat_name};
		
		my $strand = $gene_obj_ref->get_orientation();
			
			
		my $utr_seq = $gene_obj_ref->get_5prime_UTR_sequence(\$genome_seq);
		
		if ($require_5prime_UTR_flag && ! $utr_seq) { next; } # ignore those that do not have 5' UTRs.

		my ($gene_lend, $gene_rend) = sort {$a<=>$b} $gene_obj_ref->get_gene_span();

		my $promoter_region = "";

		if ($strand eq '+') {
			my $lend = $gene_lend - $upstream;
					   
			if ($lend < 1) { 
				print STDERR "warning, couldn't extract promoter length for $model, $gene_id since out of genomic range\n";
				next; # ignore entries that don't have the full required length extraction
			}

			my $upstream_seq = uc substr($genome_seq, $lend - 1, $upstream);
			my $downstream_seq = lc substr($genome_seq, $gene_lend - 1, $downstream);

			$promoter_region = $upstream_seq . $downstream_seq;
			
			unless ($genome_seq =~ /$promoter_region/i) { 
				die "Error, extracted promoter but doesn't anchor to the genome sequence from which it was extracted!";
			}

		}
		else {
			## minus strand

			my $rend = $gene_rend - $downstream;
			if ($rend < 1) {
				print STDERR "warning, couldn't extract promoter length for $model, $gene_id since out of genomic range\n";
				next;
			}
			
			my $upstream_seq = uc substr($genome_seq, $gene_rend - 1 + 1, $upstream);
			my $downstream_seq = lc substr($genome_seq, $rend, $downstream);
			
			$promoter_region = $downstream_seq . $upstream_seq;
			
			# verify
			unless ($genome_seq =~ /$promoter_region/i) { 
				die "Error, extracted promoter but doesn't anchor to the genome sequence from which it was extracted!";
			}
			
			$promoter_region = &reverse_complement($promoter_region);
		
			
		}

		unless (length $promoter_region == $upstream + $downstream) {
			print STDERR "warning, couldn't extract promoter length for $model, $gene_id since out of genomic range\n";
			next;
		}
		
		print ">$model-promoter $gene_id [$strand] upstream:$upstream downstream:$downstream\n$promoter_region\n";
		
	}
	
}



exit(0);

