#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Gene_obj_indexer;
use GFF3_utils;
use Carp;
use Nuc_translator;
use CdbTools;

$|++;

my $usage = "\n\nusage: $0 file.inx [genomeSeqFile]\n\n";

my $inx_file = $ARGV[0] or die $usage;
my $genome_seq_file = $ARGV[1] or die $usage;

unless (-s $inx_file) {
    die $usage;
}

my $gene_obj_indexer = new Gene_obj_indexer( { "use" => $inx_file } );

my $genome_seq = "";
my $prev_asmbl_id = "";


my @gene_ids = $gene_obj_indexer->get_keys();
foreach my $gene_id (sort @gene_ids) {
    my $gene_obj = $gene_obj_indexer->get_gene($gene_id);
    
    my $source = $gene_obj->{source} || ".";
    my $contig_id = $gene_obj->{asmbl_id};

    my @intron_coords = $gene_obj->get_intron_coordinates();
    my $strand = $gene_obj->get_orientation();
    my $model_id = $gene_obj->{Model_feat_name};

	if ($genome_seq_file && $contig_id ne $prev_asmbl_id) {
		$genome_seq = cdbyank_linear($contig_id, $genome_seq_file);
		$prev_asmbl_id = $contig_id;
	}


    foreach my $intron (@intron_coords) {
        my ($intron_lend, $intron_rend) = sort {$a<=>$b} @$intron;
        print "$contig_id\t$source\t$model_id\t[$strand]\t$intron_lend\t$intron_rend";

		if ($genome_seq) {
			my $intron_seq = substr($genome_seq, $intron_lend - 1, $intron_rend - $intron_lend + 1);
			if ($strand eq '-') {
				$intron_seq = &reverse_complement($intron_seq);
			}
			print "\t$intron_seq";
		}
		print "\n";
    }
}

exit(0);



