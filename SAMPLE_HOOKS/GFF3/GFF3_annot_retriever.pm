package GFF3::GFF3_annot_retriever;

use strict;
use warnings;
use Carp;
use FindBin;
use lib qw($FindBin::Bin/../PerlLib);
use base qw(PASA_UPDATES::Pasa_latest_annot_retrieval_adapter);

use GFF3_utils;
use Gene_obj_indexer;

## This sample annotation retriever is responsible for 
## parsing gene structures from a gff3 formatted file.


#### static method:
sub get_annot_retriever {
    ## the hooks is always called with a hashref for the configuration values, followed by the custom parameter.
    
    my ($config_href, $gff3_filename) = @_;

    unless (-s $gff3_filename) {
        confess "Error, cannot locate gff3 file $gff3_filename";
    }

    my $gene_obj_indexer = new Gene_obj_indexer( { create => "$gff3_filename.$$.inx"} );
    &GFF3_utils::index_GFF3_gene_objs($gff3_filename, $gene_obj_indexer);
    

    my $annot_retriever_obj = new GFF3::GFF3_annot_retriever ($gene_obj_indexer);
    
    return ($annot_retriever_obj);
    
}
   

## object methods below

####
sub new {
    my $packagename = shift;
    
    my ($gene_obj_indexer) = @_;

    my $self = $packagename->SUPER::new();
    $self->{gene_obj_indexer} = $gene_obj_indexer;

    ## get contig to gene list
    my $contig_to_genelist_aref = $self->{contig_to_genelist} = {};

    my @gene_ids = $gene_obj_indexer->get_keys();
    
    foreach my $gene_id (@gene_ids) {
        
        my $gene_obj = $gene_obj_indexer->get_gene($gene_id);

        my $contig_id = $gene_obj->{asmbl_id} or die "Error, gene_obj lacks asmbl_id value: " . $gene_obj->toString();

        my $gene_list_aref = $contig_to_genelist_aref->{$contig_id};
        unless (ref $gene_list_aref) {
            $gene_list_aref = $contig_to_genelist_aref->{$contig_id} = [];
        }

        push (@$gene_list_aref, $gene_id);
    }

    return ($self);
}

####
sub retrieve_gene_models_for_contig {
    my $self = shift;
    my ($contig_id) = @_;

    my @gene_objs;
    
    my $gene_list_aref = $self->{contig_to_genelist}->{$contig_id};
    if (ref $gene_list_aref) {
        foreach my $gene_id (@$gene_list_aref) {
            my $gene_obj = $self->{gene_obj_indexer}->get_gene($gene_id);
            unless (ref $gene_obj) {
                confess "Error, couldn't retrieve gene object for gene_id: $gene_id";
            }
            push (@gene_objs, $gene_obj);
        }
    }
    else {
		print STDERR "Warning, no genes retrieved for contig_id: $contig_id\n";
    }

    return (@gene_objs);
}



1; #EOM



