package GFF3::GFF3_annot_updater;

use strict;
use warnings;
use FindBin;
use lib qw($FindBin::Bin/../PerlLib);
use base qw(PASA_UPDATES::Pasa_update_adapter); 


## static method to retrun a Sample_annot_updater object:
sub get_updater_obj {
    my ($pasa_db_name, $custom_param) = @_;
    
    ## pasa_db_name always provided here as first parameter.
    ## custom_param can be anything (such as relational database connection params), although it's not going to be used here.
    
    my $updater_obj = new GFF3::GFF3_annot_updater($pasa_db_name);

    return ($updater_obj);
}



## object methods below

sub new {
    my $packagename = shift;

    my ($pasa_db_name) = @_;
    
    my $self = $packagename->SUPER::new($pasa_db_name);
    
    ## not going to use the custom param for anything here, 
    ## instead, just going to output updates in gff3 format.
    
    return ($self);
}


## perform single gene structure update
sub update_gene {
    my $self = shift;
    my ($annot_update_object) = @_;

    my ($updated_gene_obj) = &_print_update_to_stdout("single gene update", $annot_update_object);
    
    return ($updated_gene_obj);
    
}

## add model for an alternatively spliced variant
sub add_splicing_isoform {
    my $self = shift;
    my ($annot_update_object) = @_;
    
    my ($alt_splicing_isoform_gene_obj) = &_print_update_to_stdout("add alternative splicing isoform", $annot_update_object);
    
    return ($alt_splicing_isoform_gene_obj);
    
}

## add a new gene
sub create_new_gene {
    my $self = shift;
    my ($annot_update_object) = @_;
    
    my ($new_gene_obj) = &_print_update_to_stdout("add a new gene model", $annot_update_object);
    
    return ($new_gene_obj);

}
    
## merge several genes into a single gene model
sub merge_genes {
    my $self = shift;
    my ($annot_update_object) = @_;
    
    my ($merged_gene_obj) = &_print_update_to_stdout("merge genes", $annot_update_object);
    
    return ($merged_gene_obj);


}

## split a single gene into several models
sub split_gene {
    my $self = shift;
    my ($annot_update_object) = @_;
    
    my @split_gene_objs = &_print_update_to_stdout("split gene", $annot_update_object);
    
    return (@split_gene_objs);
    
}



###### private methods:

####
sub _print_update_to_stdout {
    my ($update_type, $annot_update_obj) = @_;
    
    my (@gene_info_structs) = $annot_update_obj->get_genes_data();
    
    my $operation_text = "# Update[$update_type] impacts the following genes: ";
    foreach my $gene_info_struct (@gene_info_structs) {
        my $gene_id = $gene_info_struct->{gene_id};
        my $model_id = $gene_info_struct->{model_id};

        $operation_text .= "(gene:$gene_id,model:$model_id),";
    }
    chop $operation_text; #remove trailing comma
    
    print "\n\n$operation_text\n";
    
    my (@gene_objs) = $annot_update_obj->get_gene_objs();
    
    foreach my $gene_obj (@gene_objs) {
        print $gene_obj->to_GFF3_format() . "\n\n";
    }
    
    
    return (@gene_objs);

}


1; #EOM

