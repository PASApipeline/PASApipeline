package PASA_UPDATES::Annot_update;

use strict;
use warnings;
use Carp;
use Gene_obj;
use Storable qw (thaw);

# generate temp identifiers for new genes outputted here.
my $gene_count = 0;
my $model_count = 0; 


sub new {
    my $packagename = shift;
    my $dbproc = shift;
    my $asmbl_id = shift;
    my $update_type = shift; # update,merge,split,alt_splice,novel

    unless ($dbproc && $asmbl_id && $update_type) {
        confess  "Error, need dbproc and asmbl_id as params";
    }

    my $self = {
        dbproc => $dbproc,
        annotdb_asmbl_id => $asmbl_id,
        update_type => $update_type,

        
        update_ids => [],
        gene_info_structs => [], 

        update_gene_objs => [],  ## if populated, used directly instead of retrieving gene_objs 
                                 ## based on the update_ids
    };

    bless ($self, $packagename);

    return ($self);
}


sub get_update_type {
    my $self = shift;
    return ($self->{update_type});
}


sub get_annotdb_asmbl_id {
    my $self = shift;
    return ($self->{annotdb_asmbl_id});
}

sub get_update_ids {
    my $self = shift;
    my @update_ids = @{$self->{update_ids}};
    return (@update_ids);
}


sub override_update_geneobjs {
    my $self = shift;
    my @update_gene_objs = @_;
    
    @{$self->{update_gene_objs}} = @update_gene_objs;
}


sub get_gene_objs {
    my $self = shift;

    my $annotdb_asmbl_id = $self->get_annotdb_asmbl_id();

    # first check for over-ride
    if (my @update_gene_objs = @{$self->{update_gene_objs}}) {
        return (@update_gene_objs);
    }

    else {

        # retrieve from the pasa database via the already stored update_ids
        
        my @update_ids = $self->get_update_ids();
        unless (@update_ids) {
            confess "Error, no update_ids stored.";
        }
        
        my @gene_objs;
        
        my $dbproc = $self->{dbproc};
        
        foreach my $update_id (@update_ids) {
            
            ## must retrieve them from the database live so as to reduce memory requirements
            my $query = "select after_gene_obj from annotation_updates where update_id = $update_id";
            my $after_gene_obj = &Mysql_connect::very_first_result_sql($dbproc, $query);
            my $new_gene_obj = thaw($after_gene_obj);
            
            die "No new gene_obj from update: $update_id " unless (ref $new_gene_obj);
            
            $new_gene_obj->{asmbl_id} = $annotdb_asmbl_id; # ensure this is set.
            $new_gene_obj->{com_name} = "updated gene (pasa update_id: $update_id)";
            
            # give temp identifiers
            $gene_count++;
            $model_count++;
            $new_gene_obj->{TU_feat_name} = "temp_gene_$gene_count";
            $new_gene_obj->{Model_feat_name} = "temp_model_$model_count";
            
            push (@gene_objs, $new_gene_obj);
        }
        
        return (@gene_objs);
    }
}

sub get_genes_data {
    my $self = shift;
    return (@{$self->{gene_info_structs}});
}


sub add_update_id {
    my $self = shift;
    my $update_id = shift;
    
    my $already_have_flag = 0;
    foreach my $already_have_update_id ($self->get_update_ids()) {
        if ($already_have_update_id == $update_id) {
            $already_have_flag = 1;
            last;
        }
    }

    unless ($already_have_flag) {
        push (@{$self->{update_ids}}, $update_id);
    }
}

sub add_gene_info {
    my $self = shift;
    my ($gene_id, $model_id) = @_;
    unless ($gene_id) {
        if (defined($model_id) && ! $gene_id) {
            $gene_id = "gene_for_" . $model_id;
        }
        else {
            confess "Error, no gene_id specified. (@_)";
        }
    }
    
    push (@{$self->{gene_info_structs}}, { gene_id =>$gene_id,
                                           model_id => $model_id } );
}


sub toString {
    my $self = shift;

    my $retText = $self->get_update_type() . ", "
        . "updateIDs: " . join (",", $self->get_update_ids()) . ",";
    
    my @gene_info_structs = $self->get_genes_data();
    foreach my $gene_info (@gene_info_structs) {
        my $gene_id = $gene_info->{gene_id};
        my $model_id = $gene_info->{model_id} || ""; # sometimes just have gene, no model id
        
        $retText .= " [gene: $gene_id, model: $model_id],";
    }

    chop $retText;

    return ($retText);
}




1; #EOM
