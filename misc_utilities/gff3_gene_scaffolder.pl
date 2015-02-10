#!/usr/bin/env perl

use strict;
use warnings;
use lib ($ENV{EUK_MODULES});
use Gene_obj;
use CdbTools;
use GFF3_utils;
use Carp;
use Nuc_translator;
use Fasta_reader;

my $usage = "usage: $0 genes.gff3 genome.fa genes.accs\n\n";

my $gff3_file = $ARGV[0] or die $usage;
my $genome_fa = $ARGV[1] or die $usage;
my $gene_accs_file = $ARGV[2] or die $usage;


my $FLANK_SIZE = 200;
my $NGAP_SIZE = 100;

main: {
    
    my %accs;
    {
        open (my $fh, $gene_accs_file) or die $!;
        while (<$fh>) {
            chomp;
            my @x = split(/\s+/);
            my $acc = $x[0];
            $accs{$acc} = 1;
        }
        close $fh;
    }
    
    my $fasta_reader = new Fasta_reader($genome_fa);
    my %genome_seqs = $fasta_reader->retrieve_all_seqs_hash();

    my $gene_obj_indexer_href = {};

    ## associate gene identifiers with contig id's.
    my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);
    
    my $fake_contig_start = $NGAP_SIZE;

    my $fake_contig = "N" x $NGAP_SIZE;

    open (my $ofh_gff3, ">$gene_accs_file.gff3") or die $!;


    foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
        my $genome_seq = $genome_seqs{$asmbl_id};
    
        my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
        foreach my $gene_id (@gene_ids) {

            unless ($accs{$gene_id}) {
                next;
            }
            
            my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
            
            print $gene_obj_ref->toString();

            my ($lend, $rend) = sort {$a<=>$b} $gene_obj_ref->get_coords();
            
            my $gene_seq = substr($genome_seq, $lend - 1 - $FLANK_SIZE, $rend - $lend + 1 + 2*$FLANK_SIZE);
            
            $fake_contig .= $gene_seq;

            my $coord_adj = -1 * $lend + $fake_contig_start + 1 + $FLANK_SIZE;
            $gene_obj_ref->adjust_gene_coordinates($coord_adj);
            
            $gene_obj_ref->{asmbl_id} = "genome";
            print $ofh_gff3 $gene_obj_ref->to_GFF3_format(); 


            $fake_contig .= "N" x $NGAP_SIZE;
            
            $fake_contig_start += length($gene_seq) + 100;

        }
    }


    close $ofh_gff3;


    open (my $ofh_scaff, ">$gene_accs_file.genome.fa") or die $!;
    $fake_contig =~ s/(\S{60})/$1\n/g;

    print $ofh_scaff ">genome\n$fake_contig";
    close $ofh_scaff;
    

    exit(0);
}


