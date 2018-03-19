package PASA_UPDATES::Pasa_latest_annot_retrieval_adapter;

use strict;
use warnings;
use Carp;

sub new {
    my $packagename = shift;
    
    my $self = { };

    bless ($self, $packagename);

    return ($self);

}


##
sub retrieve_gene_models_for_contig {
    my $self = shift;
    my ($contig_id) = @_;

    ## simply create a list of gene objects that corresond to this contig like so:
    
    # foreach gene {
    #    create a hash containing all end5 => end3 coordinates for the mRNA (ie. %mRNA)
    #    do the same for the CDS ie. %CDS    The CDS should be contained by the mRNA coords.
    #    in the case where no UTRs exist, %CDS = %mRNA
    
    #    create the gene object:
    #   
    #     my $gene_obj = new Gene_obj();
    #     $gene_obj->populate_gene_obj( \%CDS, \%mRNA );
    #
    #     add required attributes to the gene object:
    # 
    #     $gene_obj->{com_name} = "hypothetical protein"; 
    #     $gene_obj->{TU_feat_name} = "myGeneID";
    #     $gene_obj->{Model_feat_name} = "myIsoformID";
    #     $gene_obj->{asmbl_id} = $contig_id;
    #
    #     refine it (min/max coords for features are set):
    # 
    #     $gene_obj->refine_gene_object();
    # 
    #     add it to your list of gene objects:
    #
    #     push (@gene_obj_list, $gene_obj);
    #
    #  }
    #
    #  return (@gene_obj_list);
    #


    confess "Error, this class must be overridden/implemented in the subclass";
    
}

1; #EOM

