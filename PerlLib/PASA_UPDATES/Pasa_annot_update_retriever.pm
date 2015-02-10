package main;
our $SEE;

package PASA_UPDATES::Pasa_annot_update_retriever;

use strict;
use warnings;
use Ath1_cdnas;
use Mysql_connect;
use PASA_UPDATES::Annot_update;
use Carp;

sub new {
    my $packagename = shift;

    my $dbproc = shift;
    
    my $self = { _dbproc=> $dbproc,
                 _core_retrieval_query => "",
                 compare_id => undef,
             };
    
    bless ($self, $packagename);
    
    $self->_init();

    return ($self);
}


sub _init {
    my $self = shift;
    
    my $dbproc = $self->{_dbproc};
    
    my $compare_id = &Ath1_cdnas::get_max_compare_id($dbproc);
    
    $self->{compare_id} = $compare_id;

    # querying different update categories by further restricting this query to update id classes
    $self->{_core_retrieval_query} = "select distinct a.update_id, a.gene_id, a.model_id, "
        . " a.alt_splice_flag, a.is_novel_flag from annotation_updates a, status s, "
        . " status_link sl where a.compare_id = $compare_id and a.is_valid = 1 "
        . " and a.have_after = 1 and sl.compare_id = $compare_id and s.status_id = sl.status_id "
        . " and s.requires_update = 1 and sl.annot_update_id = a.update_id ";
    
}

sub get_compare_id {
    my $self = shift;
    return ($self->{compare_id});
}


sub _get_core_retrieval_query {
    my $self = shift;
    return ($self->{_core_retrieval_query});
}


sub get_single_gene_updates {
    my $self = shift;
    my ($restricted_update_ids_href,
        $restricted_status_ids_href,
        $do_not_update_ids_href)  = @_;

    unless (ref $restricted_update_ids_href eq "HASH"
            &&
            ref $restricted_status_ids_href eq "HASH"
            &&
            ref $do_not_update_ids_href eq "HASH") {
        
        confess "Error, need three hashrefs as params: @_\n";
    }
    
    my $dbproc = $self->{_dbproc};

    my %restricted_status_ids;
    my %restricted_update_ids;

    if (ref $restricted_update_ids_href) {
        %restricted_update_ids = %$restricted_update_ids_href;
    }
    if (ref $restricted_status_ids_href) {
        %restricted_status_ids = %$restricted_status_ids_href;
    }
    

    my @annot_update_objs; # populate and return

    my $query = $self->_get_core_retrieval_query() 
        . "and a.alt_splice_flag = 0 and a.is_novel_flag = 0 and s.status_id not in (29,36,40)";
    
    my @results = &Mysql_connect::do_sql_2D($dbproc, $query);
    foreach my $result (@results) {
        my ($update_id, $gene_id, $model_id, $alt_splice_flag, $novel_flag) = @$result;
    
        if ($do_not_update_ids_href->{$update_id}) {
            print "sorry, update_id: $update_id was flagged as NOT to be executed.\n";
            next;
        }
        
        
        if (%restricted_update_ids  && !$restricted_update_ids{$update_id}) {
            next; # not specified
        }
        
        
        if (%restricted_status_ids) {
            my $query = "select status_id from status_link where annot_update_id = $update_id";
            my @results = &do_sql_2D($dbproc, $query);
            my $OK = 0;
            my $status_id;
            foreach my $result (@results) {
                $status_id = $result->[0];
                if ($restricted_status_ids{$status_id}) {
                    $OK = 1;
                    print "Update status recognized, status_id = $status_id\n";
                    last;
                }
            }
            unless ($OK) {
                print "sorry, status_id($status_id) is not targeted by restricted updates.\n";
                next;
            }
        }
        
        my $asmbl_id = $self->_get_asmbl_id_via_update_id($update_id);
        print "Update to single gene model: asmbl_id($asmbl_id), update_id($update_id), gene_id($gene_id), model_id($model_id)\n";
        
        my $annot_update_obj = new PASA_UPDATES::Annot_update($dbproc, $asmbl_id, "update");
        $annot_update_obj->add_update_id($update_id);
        $annot_update_obj->add_gene_info($gene_id, $model_id);
        
        push (@annot_update_objs, $annot_update_obj);
    }

    return (@annot_update_objs);
    
}


####
sub get_split_gene_updates {
    my $self = shift;
    my ($restricted_update_ids_href,
        $do_not_update_ids_href) = @_;

    unless (ref $restricted_update_ids_href eq "HASH"
            &&
            ref $do_not_update_ids_href eq "HASH") {
        confess "Error, params require two hashrefs";
    }
    
    my $dbproc = $self->{_dbproc};
    
    my $query = $self->_get_core_retrieval_query() 
        . " and a.alt_splice_flag = 0 and a.is_novel_flag = 0 and s.status_id = 40 order by a.gene_id";
    

    ## group split gene objects by the original model identifier:
    my %model_id_to_new_model_info;
    my @results = &Mysql_connect::do_sql_2D($dbproc, $query);
    foreach my $result (@results) {
        my ($update_id, $gene_id, $model_id, $alt_splice_flag, $is_novel_flag) = @$result;
        
        if ($do_not_update_ids_href->{$update_id}) {
            print "Sorry, splitting provided by update_id: $update_id has been indicated NOT to occur\n";
            next;
        }

        if (%$restricted_update_ids_href && ! $restricted_update_ids_href->{$update_id}) {
            print "Sorry, targeting only specific update_ids and update_id: $update_id was not flagged for execution.\n";
            next;
        }
        

        my $update_struct = { update_id => $update_id,
                              gene_id => $gene_id,
                              model_id => $model_id};
        my $list_aref = $model_id_to_new_model_info{$model_id};
        unless (ref $list_aref) {
            $list_aref = $model_id_to_new_model_info{$model_id} = [];
        }
        push (@$list_aref, $update_struct);
    }
    
    
    ## Check each model, make sure more than one entry exists for each.
    foreach my $model_id (keys %model_id_to_new_model_info) {
        my $list_aref = $model_id_to_new_model_info{$model_id};
        my $num_genes = scalar (@$list_aref);
        if ($num_genes < 2) {
            confess "Error, $model_id model is supposedly being split but has less than two new gene models!\n";
        }
    }
    
    
    
    ## Prepare annot update objs:

    my @annot_update_objs;

    
    foreach my $model_id (keys %model_id_to_new_model_info) {
        print "Gene Splitting: $model_id\n" if $main::SEE;
        my $list_aref = $model_id_to_new_model_info{$model_id};
        ## perform update on first gene, then add others as new genes:
        my $first_entry = shift @$list_aref;
        my ($update_id, $gene_id, $model_id) = ($first_entry->{update_id}, $first_entry->{gene_id}, $first_entry->{model_id});
        my $gene_model_being_split = "$gene_id,$model_id";
        print "\tupdating model as part of gene split: update_id($update_id), gene_id($gene_id), model_id($model_id)\n" if $main::SEE;
        
        my $asmbl_id = $self->_get_asmbl_id_via_update_id($update_id);
        my $annot_update_obj = new PASA_UPDATES::Annot_update($dbproc, $asmbl_id, "split");
        $annot_update_obj->add_update_id($update_id);
        $annot_update_obj->add_gene_info($gene_id, $model_id);
 
        ## get the update_ids that provide the structures for the other split products:
        while (my $entry = shift @$list_aref) {
            my ($update_id) = $entry->{update_id};
            $annot_update_obj->add_update_id($update_id);
        }
        
        push (@annot_update_objs, $annot_update_obj);
    }


    return (@annot_update_objs);
}
        


sub get_merge_gene_updates {
    my $self = shift;
    
    my ($restricted_update_ids_href,
        $do_not_update_ids_href) = @_;
    
    unless (ref $restricted_update_ids_href eq "HASH"
            &&
            ref $do_not_update_ids_href eq "HASH") {
        confess "Error, params require two hashrefs";
    }
    
    

    my $dbproc = $self->{_dbproc};
    
    my $query = $self->_get_core_retrieval_query();
    $query .= " and a.alt_splice_flag = 0 and a.is_novel_flag = 0 and s.status_id in (29,36)";
    
    my @annot_update_objs;
    

    my @results = &Mysql_connect::do_sql_2D($dbproc, $query);
    foreach  my $result (@results) {
        my ($update_id, $gene_id, $model_id, $alt_splice_flag, $is_novel_flag) = @$result;
        print "Gene Merge. update_id($update_id), gene_id($gene_id), model_id($model_id)\n";
        
        if ($do_not_update_ids_href->{$update_id}) {
            print "Sorry, merging operation provided by update_id: $update_id has been indicated NOT to occur\n";
            next;
        }
        
        if (%$restricted_update_ids_href && ! $restricted_update_ids_href->{$update_id}) {
            print "Sorry, targeting only specific update_ids and update_id: $update_id was not flagged for execution.\n";
            next;
        }
        

        my @gene_ids = split (/\s+/, $gene_id);
                
        unless (scalar (@gene_ids) > 1) {
            confess "Error, update_id: $update_id should be a merging event but lacks multiple genes to merge: $gene_id";
        }
        
        my $asmbl_id =  $self->_get_asmbl_id_via_update_id($update_id);

        my $annot_update_obj = new PASA_UPDATES::Annot_update($dbproc, $asmbl_id, "merge");
        
        $annot_update_obj->add_update_id($update_id);
        foreach my $gene_id (@gene_ids) {
            
            $annot_update_obj->add_gene_info($gene_id, undef); #don't provide models here.  The gene and model list are unlikely to corrspond in their ordering, unfortunately.
        }
        
        push (@annot_update_objs, $annot_update_obj);

    }

    return (@annot_update_objs);
    
}


sub get_novel_gene_updates {
    my $self = shift;
    my ($restricted_update_ids_href,
        $do_not_update_ids_href) = @_;
    
    unless (ref $restricted_update_ids_href eq "HASH"
            &&
            ref $do_not_update_ids_href eq "HASH") {
        confess "Error, params require two hashrefs";
    }
    
    my $dbproc = $self->{_dbproc};
    
    my $query = $self->_get_core_retrieval_query() 
        . " and is_novel_flag = 1 and alt_splice_flag = 0 and is_valid = 1 ";
    
    
    my @annot_update_objs;


    my @results = &Mysql_connect::do_sql_2D($dbproc, $query);
    foreach my $result (@results) {
        my ($update_id, $gene_id, $model_id, $alt_splice_flag, $novel_flag) = @$result;
        
        if ($do_not_update_ids_href->{$update_id}) {
            print "Sorry, novel gene provided by update_id: $update_id has been indicated NOT to occur\n";
            next;
        }

        if (%$restricted_update_ids_href && ! $restricted_update_ids_href->{$update_id}) {
            print "Sorry, targeting only specific update_ids and update_id: $update_id was not flagged for execution.\n";
            next;
        }
        
        


        my $asmbl_id =  $self->_get_asmbl_id_via_update_id($update_id);

        my $annot_update_obj = new PASA_UPDATES::Annot_update($dbproc, $asmbl_id, "novel");

        $annot_update_obj->add_update_id($update_id);
        $annot_update_obj->add_gene_info($gene_id, $model_id);
        
        push (@annot_update_objs, $annot_update_obj);
    }
        
    return (@annot_update_objs);
}


sub get_alt_splice_isoform_updates {
    my $self = shift;
    my ($restricted_update_ids_href,
        $do_not_update_ids_href) = @_;
    
    unless (ref $restricted_update_ids_href eq "HASH"
            &&
            ref $do_not_update_ids_href eq "HASH") {
        confess "Error, params require two hashrefs";
    }
    
    my $dbproc = $self->{_dbproc};
    
    my $query = $self->_get_core_retrieval_query() 
        . " and is_valid = 1 and alt_splice_flag = 1";
    
    
    my @annot_update_objs;


    my @results = &Mysql_connect::do_sql_2D($dbproc, $query);
    foreach my $result (@results) {
        my ($update_id, $gene_id, $model_id, $alt_splice_flag, $novel_flag) = @$result;
        
        if ($do_not_update_ids_href->{$update_id}) {
            print "Sorry, altsplice isoform provided by update_id: $update_id has been indicated NOT to occur\n";
            next;
        }
        
        if (%$restricted_update_ids_href && ! $restricted_update_ids_href->{$update_id}) {
            print "Sorry, targeting only specific update_ids and update_id: $update_id was not flagged for execution.\n";
            next;
        }
        
        print "-alt_splice isoform added: update_id($update_id), gene_id($gene_id), model_id($model_id)\n";
        
        my $asmbl_id =  $self->_get_asmbl_id_via_update_id($update_id);

        my $annot_update_obj = new PASA_UPDATES::Annot_update($dbproc, $asmbl_id, "alt_splice");
        
        $annot_update_obj->add_update_id($update_id);
        $annot_update_obj->add_gene_info($gene_id, $model_id);
        
        push (@annot_update_objs, $annot_update_obj);
    }

    return (@annot_update_objs);
}


sub _get_asmbl_id_via_update_id {
    my $self = shift;
    my $update_id = shift;
    
    my $dbproc = $self->{_dbproc};
    
    my $query = "select cdna_acc from status_link where annot_update_id = $update_id";
    my $cdna_acc = &Mysql_connect::very_first_result_sql($dbproc, $query);
    
    $query = "select c.annotdb_asmbl_id from clusters c, cluster_link cl where cl.cdna_acc = ? and cl.cluster_id = c.cluster_id";
    my $asmbl_id = &Mysql_connect::very_first_result_sql($dbproc, $query, $cdna_acc);
    return ($asmbl_id);
}



1; #EOM

