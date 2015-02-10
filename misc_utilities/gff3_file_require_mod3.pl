#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Fasta_reader;
use GFF3_utils;
use Carp;
use Nuc_translator;
use Longest_orf;
use Data::Dumper;

my $usage = "\n\nusage: $0 gff3_file genome_db [genetic_code=universal|Tetrahymena|Candida]\n\n";

my $gff3_file = $ARGV[0] or die $usage;
my $fasta_db = $ARGV[1] or die $usage;
my $genetic_code = $ARGV[2];


if ($genetic_code) {
    &Nuc_translator::use_specified_genetic_code($genetic_code);
}

my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

my $fasta_reader = new Fasta_reader($fasta_db);
my %genome = $fasta_reader->retrieve_all_seqs_hash();

my $TOTAL_MODELS_ADJUSTED = 0;

foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
    my $genome_seq = $genome{$asmbl_id} or die "Error, cannot locate sequence for genome entry: $asmbl_id";
    
    my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
    foreach my $gene_id (@gene_ids) {
        my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
				
        foreach my $isoform ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
 
		
            my $cds_length = $isoform->get_CDS_length();
            if ($cds_length % 3 == 0) {
                next;
            }

			my $orientation = $isoform->get_orientation();
			my ($model_lend, $model_rend) = sort {$a<=>$b} $isoform->get_model_span();
			my ($gene_lend, $gene_rend) = sort {$a<=>$b} $isoform->get_gene_span();
			
            my $isoform_id = $isoform->{Model_feat_name};
            
            my $cds_seq = $isoform->create_CDS_sequence(\$genome_seq);
            
            #print "$isoform_id\t$cds_seq\n";
            
            my $longest_orf = new Longest_orf();
            $longest_orf->allow_partials();
            $longest_orf->forward_strand_only();
            
            my $orf = $longest_orf->get_longest_orf($cds_seq);
            
            #print Dumper($orf);
            
            my $orf_start = $orf->{start};
            my $orf_end = $orf->{stop};
            if ($orf_end > $cds_length) {
                $orf_end -= 3;
            }
            
            my $coding_len = length($cds_seq);
            
            #print "$cds_length vs. coding $coding_len\n";
            
            my $adjusted_terminus = 0;
            if ($orf_start != 1 && $orf_start < 3) {
                &adjust_orf_start($isoform, $orf_start-1);
                $adjusted_terminus++;
            }
            if ($orf_end != $cds_length && $cds_length - $orf_end < 3) {
                &adjust_orf_end($isoform, $cds_length - $orf_end);
                $adjusted_terminus++;
            }
            
            if ($adjusted_terminus) {
                $TOTAL_MODELS_ADJUSTED++;
            }
            

        }
        
        print $gene_obj_ref->to_GFF3_format() . "\n";
        

    }
    

}


print "\n\n# Total models adjusted: $TOTAL_MODELS_ADJUSTED\n\n";

exit(0);

####
sub adjust_orf_start {
    my ($isoform, $delta) = @_;

    my $model_name = $isoform->{Model_feat_name};

    my $orientation = $isoform->get_orientation();

    print STDERR "$model_name\tAdjust start by $delta\n";

    my @exons = $isoform->get_exons();
    
    while (@exons) {
        my $exon = shift @exons;
        if (my $cds = $exon->get_CDS_exon_obj()) {
            if ($orientation eq '+') {
                $cds->{end5} += $delta;
                last;
            }
            else {
                $cds->{end5} -= $delta;
                last;
            }
        }
    }
    
    $isoform->refine_gene_object();
    

    

}


sub adjust_orf_end {
    my ($isoform, $delta) = @_;

    my $model_name = $isoform->{Model_feat_name};

    my $orientation = $isoform->get_orientation();

    print STDERR "#$model_name\tAdjust stop by $delta\n";
    
    my @exons = reverse $isoform->get_exons();
    
    while (@exons) {
        my $exon = shift @exons;
        if (my $cds = $exon->get_CDS_exon_obj()) {
            if ($orientation eq '+') {
                $cds->{end3} -= $delta;
                last;
            }
            else {
                $cds->{end3} += $delta;
                last;
            }
        }
    }
    
    $isoform->refine_gene_object();
    

    return;
}

