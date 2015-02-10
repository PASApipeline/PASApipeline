package PASA_UPDATES::Pasa_update_adapter;

use strict;
use warnings;
use Carp;


sub new {
    my $packagename = shift;
    my ($pasa_db_name) = @_;
    
    my $self = {
        pasa_db_name => $pasa_db_name,
        debug => 0,
        
    };

    bless ($self, $packagename);

    return ($self);
}


## perform single gene structure update
sub update_gene {
    my $self = shift;
    my ($annot_update_object) = @_;


    confess "Must override this method in subclass";
    
}

## add model for an alternatively spliced variant
sub add_splicing_isoform {
    my $self = shift;
    my ($annot_update_object) = @_;

    confess "Must override this method in subclass";
    
}

## add a new gene
sub create_new_gene {
    my $self = shift;
    my ($annot_update_object) = @_;
    
    confess "Must override this method in subclass";

}
    
## merge several genes into a single gene model
sub merge_genes {
    my $self = shift;
    my ($annot_update_object) = @_;

    confess "Must override this method in subclass";
}

## split a single gene into several models
sub split_gene {
    my $self = shift;
    my ($annot_update_object) = @_;

    confess "Must override this method in subclass";
}


1; #EOM


    
