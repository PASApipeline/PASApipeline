#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Gene_obj_indexer;
use CdbTools;
use GFF3_utils;
use Carp;
use Fasta_reader;
use Data::Dumper;

my $usage = "\n\nusage: $0 gff3_file genome_db [gene|exon|CDS]\n\n";

my $gff3_file = $ARGV[0] or die $usage;
my $genome_db = $ARGV[1] or die $usage;
my $seq_type = $ARGV[2] || "gene";


unless ($seq_type =~ /^(gene|exon|CDS)$/) {
    die "Error, don't understand sequence type [$seq_type]\n\n$usage";
}

my $index_file = "$gff3_file.inx";

my $gene_obj_indexer = undef;

if (-s $index_file) {
    ## try to use it
    $gene_obj_indexer = new Gene_obj_indexer( { "use" => $index_file } );
    my @gene_ids = $gene_obj_indexer->get_keys();
    unless (@gene_ids) {
        $gene_obj_indexer = undef; # didn't work, must create a new index file
        print STDERR "Even though $index_file exists, couldn't use it.  Going to have to regenerate it now.\n";
    }
}

unless ($gene_obj_indexer) {
    
    $gene_obj_indexer = new Gene_obj_indexer( { "create" => $index_file } );
    
    ## associate gene identifiers with contig id's.
    &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer);
}

## associate all gene_ids with contigs
my @all_gene_ids = $gene_obj_indexer->get_keys();
my %contig_to_gene_list;
foreach my $gene_id (@all_gene_ids) {
    my $gene_obj = $gene_obj_indexer->get_gene($gene_id);
    
    my $contig = $gene_obj->{asmbl_id} 
    or croak "Error, can't find contig id attached to gene_obj($gene_id) as asmbl_id val\n" 
        . $gene_obj->toString();
    
    
    my $gene_list_aref = $contig_to_gene_list{$contig};
    unless (ref $gene_list_aref) {
        $gene_list_aref = $contig_to_gene_list{$contig} = [];
    }

    push (@$gene_list_aref, $gene_id);

}


## mask genes
my $fasta_reader = new Fasta_reader($genome_db);
while (my $seq_obj = $fasta_reader->next()) {
	
	my $accession = $seq_obj->get_accession();
	my $header= $seq_obj->get_header();
	my $sequence = lc $seq_obj->get_sequence();
	
	my @seq_chars = split (//, $sequence);

	my $gene_list_aref = $contig_to_gene_list{$accession};
	
	foreach my $gene_id (@$gene_list_aref) {
		print STDERR "-masking $gene_id on $accession\n";
		my $gene_obj = $gene_obj_indexer->get_gene($gene_id);

		my @coords = &retrieve_relevant_coords($gene_obj);
		
		# print STDERR "-masking coords: " . Dumper (\@coords);

		foreach my $coordset (@coords) {
			my ($lend, $rend) = sort {$a<=>$b} @$coordset;
			for (my $i = $lend; $i <= $rend; $i++) {
				$seq_chars[$i-1] = 'N';
			}
		}
	}

	# output fasta format:
	$sequence = join ("", @seq_chars);
	$sequence =~ s/(\w{60})/$1\n/g;
	chomp $sequence;

	print ">$header ($seq_type masked)\n$sequence\n";
	
	
}


exit(0);



####
sub retrieve_relevant_coords {
	my ($gene_obj) = @_;

	my @coords;

	if ($seq_type eq 'gene') {
		my ($end5, $end3) = $gene_obj->get_coords();
		push (@coords, [$end5, $end3]);
	}
	else {
		my @exons = $gene_obj->get_exons();
		foreach my $exon (@exons) {
			my ($end5, $end3) = $exon->get_coords();
			push (@coords, [$end5, $end3]) if ($seq_type eq 'exon');
			if ($seq_type eq 'CDS' && (my $cds = $exon->get_CDS_exon_obj()) ) {
				my ($end5, $end3) = $cds->get_coords();
				push (@coords, [$end5, $end3]);
			}
		}
	}

	foreach my $isoform ($gene_obj->get_additional_isoforms()) {
		my @other_coords = &retrieve_relevant_coords($isoform);
		push (@coords, @other_coords);
	}
	
	return (@coords);
}
