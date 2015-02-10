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

my $usage = "\n\nusage: $0 gff3_file genome_db core_outputfilename\n\n";

my $gff3_file = $ARGV[0] or die $usage;
my $fasta_db = $ARGV[1] or die $usage;
my $core_outputfile = $ARGV[2] or die $usage;

my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);


## create various output files:

open (my $exons_fh, ">$core_outputfile.exons") or die "Error, cannot write to file $core_outputfile.exons";
open (my $introns_fh, ">$core_outputfile.introns") or die "Error, cannot write to file $core_outputfile.introns";


my $num_exons = 0;
my $num_genes = 0;


foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
    my $genome_seq = cdbyank_linear($asmbl_id, $fasta_db);
    
    my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
    foreach my $gene_id (@gene_ids) {
        		
		$num_genes++;
		
		my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
		
                
		my $TU_feat_name = $gene_obj_ref->{TU_feat_name};
		my $model_feat_name = $gene_obj_ref->{Model_feat_name};
		
		print STDERR "-processing $TU_feat_name\t$model_feat_name\n";
		
		my $orientation = $gene_obj_ref->get_orientation();
		
		
		my @exons = $gene_obj_ref->get_exons();
		
		foreach my $exon (@exons) {
			my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
			my $exon_length = $rend - $lend +1;
			my $exon_seq = substr($genome_seq, $lend - 1, $rend - $lend + 1);
			if ($orientation eq '-') {
				$exon_seq = &reverse_complement($exon_seq);
			}

			my $percent_gc = &calc_percent_gc($exon_seq);
			
			print $exons_fh "$asmbl_id\t$TU_feat_name\t$model_feat_name\t$lend\t$rend\t$orientation\t$exon_length\t$percent_gc\t$exon_seq\n";
			
			$num_exons++;
		}

		my @introns = $gene_obj_ref->get_intron_coordinates();
		
		foreach my $intron (@introns) {
			my ($lend, $rend) = sort {$a<=>$b} @$intron;
			my $intron_length = $rend - $lend + 1;
			my $intron_seq = substr($genome_seq, $lend - 1, $rend - $lend + 1);
			if ($orientation eq '-') {
				$intron_seq = &reverse_complement($intron_seq);
			}
			
			my $percent_gc = &calc_percent_gc($intron_seq);
			
			print $introns_fh "$asmbl_id\t$TU_feat_name\t$model_feat_name\t$lend\t$rend\t$orientation\t$intron_length\t$percent_gc\t$intron_seq\n";
		}
		
    }
}

my $num_exons_per_gene = $num_exons / $num_genes;
print "\n\nNumber of exons per gene: $num_exons_per_gene\n";


exit(0);


####
sub calc_percent_gc {
	my ($seq) = @_;

	$seq =~ s/n//ig;
	
	my $length = length($seq);
	
	unless ($length > 0) {
		print STDERR "Error, trying to compute %GC with sequence of zero length!\n";
		return(0);
	}
	
	my $gc_count = 0;
	while ($seq =~ /[gc]/ig) {
		$gc_count++;
	}
	
	my $percent_gc = sprintf("%.2f", $gc_count / $length * 100);
	
	return($percent_gc);
}


