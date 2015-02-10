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
use Nuc_translator;

my $usage = "\n\nusage: $0 gff3_file genome_db flank_region_size\n\n";

my $gff3_file = $ARGV[0] or die $usage;
my $fasta_db = $ARGV[1] or die $usage;
my $flank_region_size = $ARGV[2];

unless (defined $flank_region_size && $flank_region_size =~ /^\d+$/) {
	die $usage;
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



foreach my $asmbl_id (sort keys %contig_to_gene_list) {
    
    my $genome_seq = cdbyank_linear($asmbl_id, $fasta_db);
    
	my $genome_seq_length = length($genome_seq);

    my @gene_ids = @{$contig_to_gene_list{$asmbl_id}};
    
    foreach my $gene_id (@gene_ids) {
        my $gene_obj_ref = $gene_obj_indexer->get_gene($gene_id);
		
        my $orientation = $gene_obj_ref->get_orientation();
		my ($lend, $rend) = sort {$a<=>$b} $gene_obj_ref->get_coords();

		$lend -= $flank_region_size;
		$lend = 1 if $lend < 1;
		
		$rend += $flank_region_size;
		$rend = $genome_seq_length if $rend > $genome_seq_length;

		my $region_seq = substr($genome_seq, $lend - 1, $rend - $lend + 1);

		if ($orientation eq '-') {
			$region_seq = reverse_complement($region_seq);
		}

		$region_seq =~ s/(\S{60})/$1\n/g;
		chomp $region_seq;

		print ">$gene_id ($lend-$rend; plus flank=$flank_region_size)\n$region_seq\n";
	}
}


exit(0);

