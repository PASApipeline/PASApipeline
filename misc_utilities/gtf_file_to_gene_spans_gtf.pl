#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Gene_obj_indexer;
use GTF_utils;
use Carp;
use Data::Dumper;

my $usage = "\n\nusage: $0 gtf_file\n\n";

my $gtf_file = $ARGV[0] or die $usage;


my $index_file = "$gtf_file.inx";

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
    &GTF_utils::index_GTF_gene_objs($gtf_file, $gene_obj_indexer);
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
    
    my @gene_ids = @{$contig_to_gene_list{$asmbl_id}};
    
    foreach my $gene_id (@gene_ids) {
        my $gene_obj_ref = $gene_obj_indexer->get_gene($gene_id);
        
        my ($gene_lend, $gene_rend) = sort {$a<=>$b} $gene_obj_ref->get_coords();
        my $strand = $gene_obj_ref->get_orientation();
        

        my $gene_id = $gene_obj_ref->{TU_feat_name};
        my $trans_id = $gene_obj_ref->{Model_feat_name};
        my $contig_id = $gene_obj_ref->{asmbl_id};

        
        print join("\t", $contig_id, "ref_annot",
                   "exon",
                   $gene_lend, $gene_rend,
                   ".",
                   $strand,
                   ".",
                   "gene_id \"$gene_id\"; transcript_id \"$trans_id\"") . "\n";
        

    }
}



exit(0);


