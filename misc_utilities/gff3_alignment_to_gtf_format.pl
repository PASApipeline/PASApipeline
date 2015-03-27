#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use GFF3_alignment_utils;
use Fasta_retriever;

my $usage = "usage: $0 trans_align.gff3\n\n";

my $trans_align_gff3 = $ARGV[0] or die $usage;


main: {

    my %cdna_alignments;
    
    print STDERR "-parsing $trans_align_gff3\n";
    my %scaff_to_cdna_ids = &GFF3_alignment_utils::index_alignment_objs($trans_align_gff3, \%cdna_alignments);


    print STDERR "-outputting GTF format\n";
    foreach my $scaff (keys %scaff_to_cdna_ids) {

        foreach my $cdna_id (@{$scaff_to_cdna_ids{$scaff}}) {

            my $cdna_obj = $cdna_alignments{$cdna_id};

            print $cdna_obj->to_GTF_format( gene_id => $cdna_obj->{gene_id} ) . "\n";
            
        }
    }


    exit(0);
}


