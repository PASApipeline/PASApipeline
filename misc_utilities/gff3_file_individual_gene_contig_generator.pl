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

my $usage = "\n\nusage: $0 gff3_file genome_db [output_prefix=gene_contigs]\n\n";

my $gff3_file = $ARGV[0] or die $usage;
my $fasta_db = $ARGV[1] or die $usage;
my $out_prefix = $ARGV[2] || "gene_contigs";


my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

my $fake_contig = "";

my $pad = 300;



open (my $ofh_seq, ">$out_prefix.genome") or die $!;
open (my $ofh_genes, ">$out_prefix.gff3") or die $!;


foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
	my $counter = 0;

    my $genome_seq = cdbyank_linear($asmbl_id, $fasta_db);
    
	unless ($genome_seq) {
		die "Error, no genome_seq for $asmbl_id";
	}


    my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
    foreach my $gene_id (@gene_ids) {
        my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
		
        my %exon_coords;
		
		my ($gene_lend, $gene_rend) = sort {$a<=>$b} $gene_obj_ref->get_coords();
		my $adjust = $gene_lend - $pad;
		if ($adjust < 0) {
			$adjust = 0;
		}
		
		#print "Before: " . $gene_obj_ref->to_GFF3_format();

		#print "\n\nAdjust: $adjust\n";
		
		$gene_obj_ref->adjust_gene_coordinates(-1 * $adjust);
							
		my $subseq = substr($genome_seq, $adjust, $gene_rend - $gene_lend + 2*$pad);
		
		$counter++;
		print $ofh_seq ">$asmbl_id.$counter\n$subseq\n";
		
		$gene_obj_ref->{asmbl_id} = "$asmbl_id.$counter";
		print $ofh_genes $gene_obj_ref->to_GFF3_format() . "\n";
	}
}
close $ofh_seq;
close $ofh_genes;

print STDERR "-done.  See files $out_prefix.gff3 and $out_prefix.genome\n\n";


exit(0);



