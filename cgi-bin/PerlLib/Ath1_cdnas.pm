#!/usr/local/bin/perl

package main;
our $SEE;

## Simple API to access the Mysql Ath1_cdnas database.

package Ath1_cdnas;
use strict;
use Mysql_connect;
use CDNA::CDNA_alignment;
use CDNA::Alignment_segment;
use Gene_obj;
use Storable qw(thaw);
use Carp;

sub get_cluster_ids_via_annotdb_asmbl_id {
    my ($dbproc, $asmbl_id) = @_;
    my $query = qq {
        select distinct c.cluster_id 
            from clusters c, align_link al, cdna_info ci
            where c.annotdb_asmbl_id = ?
            and c.cluster_id = al.cluster_id
            and ci.id = al.cdna_info_id
            and ci.is_assembly = 0
            and al.validate = 1
        };
    
    my @cluster_id_refs = &Mysql_connect::do_sql_2D ($dbproc, $query, $asmbl_id);
    my @cluster_ids;
    foreach my $cluster_ref (@cluster_id_refs) {
        my $cluster_id = $cluster_ref->[0];
        push (@cluster_ids, $cluster_id);
    }
    return (@cluster_ids);
}



####
sub create_alignment_obj {
    my $dbproc = shift;
    my $align_id = shift;
    my $seq_ref = shift;
    
    
    my $query = "select c.annotdb_asmbl_id from clusters c, align_link al where al.align_id = $align_id and al.cluster_id = c.cluster_id";
    my $genome_acc = &Mysql_connect::very_first_result_sql($dbproc, $query);


    ## get the alignment segment coordinates
    $query = "select ci.id, ci.length, ci.cdna_acc, ci.is_fli, "
        . " ci.header, a.lend, a.rend, a.mlend, a.mrend, a.orient, "
        . " al.aligned_orient, al.spliced_orient, al.validate, al.prog, al.align_acc, a.per_id "
        . " from alignment a, cdna_info ci, align_link al "
        . " where a.align_id = $align_id and a.align_id = al.align_id "
        . " and al.cdna_info_id = ci.id ";
    
    my @results = &Mysql_connect::do_sql_2D ($dbproc, $query);
    my @alignment_segments;
    my $cdna_length;
    my $Cdna_acc,
    my $Align_acc;
    my $is_fli;
    my $spliced_orientation;
    my $aligned_orientation;
    my $title;
    my $validation_status = 0;
    my $Cdna_id = undef;
    my $Prog;
    foreach my $result (@results) {
        my ($cdna_id, $length, $acc, $is_fli_status, $header, $lend, $rend, $mlend, $mrend, $orient, $aligned_orient, $spliced_orient, $validate, $prog, $align_acc, $per_id) = @$result;
        $cdna_length = $length;
        $Cdna_acc = $acc;
        $Align_acc = $align_acc,
        $is_fli = $is_fli_status;
        $title = $header;
        $spliced_orientation = $spliced_orient;
        $aligned_orientation = $aligned_orient;
        $Cdna_id = $cdna_id;
        $Prog = $prog;

        my $seg_length = abs($mrend - $mlend) + 1;
        
        if ($validate) {
            $validation_status = $validate;
        }
        
        if ($seg_length > 1 && $aligned_orient ne $orient) {
            warn "Error, $acc has alignment.orient = $orient but cdna_link aligned_orient = $aligned_orient ";
        }
        
        if ($orient eq '-') {
            ($lend, $rend) = ($rend, $lend); #swap coordinates.
        }
        my $segment = new CDNA::Alignment_segment ($lend, $rend, $mlend, $mrend, $per_id);
        push (@alignment_segments, $segment);
    }

    unless (@alignment_segments) {
        confess "Error, no alignment segments extracted for align_id: $align_id\n";
    }
    
    my $alignment = new CDNA::CDNA_alignment ($cdna_length, \@alignment_segments,  $seq_ref);
    $alignment->set_align_acc($Align_acc);
    $alignment->set_cdna_acc($Cdna_acc);
    $alignment->set_fli_status($is_fli);
    
    $alignment->set_cdna_id($Cdna_id);
    $alignment->set_align_id($align_id);
    $alignment->{genome_acc} = $genome_acc;
    $alignment->{prog} = $Prog;
    
    if ($title) {
        $alignment->set_title($title);
    }
    
    my $orientation = $alignment->get_orientation();
    unless ($orientation =~ /[\+\-]/) {
        $alignment->set_orientation('+'); ## set a default
    }
    
    my $computed_spliced_orient = $alignment->get_spliced_orientation();
    

    if (! $seq_ref) { # if no sequence, must rely on the database
        
        # old way
        #if ($spliced_orientation =~ /^[+-]$/) { ## set in database, set here explicitly, otherwise, relying on sequence analysis done during obj instantiation
        
        # now, if got seq_ref, use it, so the consensus splice sites will be populated properly.
        # otherwise, rely on the database settings.

        unless ($spliced_orientation) {
            confess "Error, no sequence as input parameter and spliced validation not avail in database.";
        }
        
        if ($spliced_orientation ne '?') {
            $alignment->set_spliced_orientation($spliced_orientation);
            $alignment->force_spliced_validation($spliced_orientation);
        }
    }
    
    else {
        ## have sequence, but still rely on database for spliced orient if it's different and non-ambiguous in the db.
        
        if ($validation_status) { #rely on spliced orient stored in db.
            if ($computed_spliced_orient eq '?' && $spliced_orientation ne '?') {
                $alignment->set_spliced_orientation($spliced_orientation);
            }
            
            if ($computed_spliced_orient ne '?' && $spliced_orientation ne $computed_spliced_orient) {
                ## huge problem. Something wrong with calculated spliced orient?
                confess "Error, spliced orient in db ($spliced_orientation) differs from calculated spliced orient ($computed_spliced_orient) ";
            }
        }
        else { 
            ## invalid alignment, but got sequence.  Rely on sequence analysis, computed spliced alignment
            ;
        }
    }
    
    return ($alignment);
}




####
sub get_alignment_obj_via_align_acc {
    my ($dbproc, $align_acc, $seqref) = @_; #seqref is optional
    
    my $align_id = &get_align_id_via_align_acc($dbproc, $align_acc);
    
    my $alignment = &create_alignment_obj($dbproc, $align_id, $seqref);
    return ($alignment);
}



####
sub get_align_id_via_align_acc {
    my ($dbproc, $align_acc) = @_;
    my $query = "select align_id from align_link where align_acc = ? ";
    my $align_id = &very_first_result_sql($dbproc, $query, $align_acc);
    unless (defined $align_id) {
        confess "Error, no align_id returned for align_acc: $align_acc";
    }

    return($align_id);
}




####
sub get_validating_align_ids_via_acc {
    my $dbproc = shift;
    my $acc = shift;
    my $query = "select align_id from align_link al, cdna_info ci "
        . " where ci.cdna_acc = ? and ci.id = al.cdna_info_id and al.validate = 1";
    my @results = &Mysql_connect::do_sql_2D($dbproc, $query, $acc);
    
    my @align_ids;
    foreach my $result (@results) {
        my $align_id = $result->[0];
        push (@align_ids, $align_id);
    }

    return(@align_ids);
}

####
sub get_cdna_info_via_align_id {
    my $dbproc = shift;
    my $align_id = shift;
    my $query = "select al.align_id, al.align_acc, al.prog, al.validate, "
        . " al.spliced_orient, al.num_segments, ci.cdna_acc, ci.is_fli, ci.header, ci.length, "
        . " al.avg_per_id, al.percent_aligned, al.alignment "
        . " from align_link al, cdna_info ci "
        . " where al.align_id = ? and al.cdna_info_id = ci.id";
    
    my $result = &Mysql_connect::first_result_sql($dbproc, $query, $align_id);
    my ($qalign_id, $align_acc, $prog, $validate, $spliced_orient, $num_segments, $cdna_acc, $is_fli, $header, $length, $avg_per_id, $percent_aligned, $alignment) = @$result;
    my $s = {align_id => $align_id,
             cdna_acc => $cdna_acc,
             align_acc => $align_acc,
             prog => $prog,
             validate => $validate,
             spliced_orient => $spliced_orient,
             num_segments => $num_segments,
             is_fli => $is_fli,
             header=>$header,
             length=>$length,
             avg_per_id=>$avg_per_id,
             percent_aligned=>$percent_aligned,
             alignment=>$alignment};
    return ($s);
}

####
sub get_subcluster_id_for_cdna_acc {
	my ($dbproc, $cdna_acc) = @_;

	my $query = "select subcluster_id from subcluster_link where cdna_acc = ?";
	my $result = &very_first_result_sql($dbproc, $query, $cdna_acc);
	return($result);
}



####
sub get_annotations_linked_to_assembly {
    my $dbproc = shift;
    my $acc = shift;
    my $compare_id = shift;
    my $query = "select gene_id, model_id from annotation_link where cdna_acc = ? and compare_id = ?";
    my @results = &do_sql_2D ($dbproc, $query, $acc, $compare_id);
    my @annots;
    foreach my $result (@results) {
        my $annot = $result->[0];
        push (@annots, [@$result]);
    }
    return (@annots);
}


####
sub get_status_info {
    my $dbproc = shift;
    my $acc = shift;
    my $compare_id = shift;
    my $query = "select sl.annot_update_id, sl.status_link_id, sl.cdna_acc, sl.status_id, sl.is_chromo, sl.comment, sl.curated_exception, sl.curator_comment, s.status_descr from status_link sl, status s where cdna_acc = ? and sl.status_id = s.status_id and sl.compare_id = ?";
    my @results = &do_sql_2D ($dbproc, $query, $acc, $compare_id);
    my @x;
    foreach my $result (@results) {
        my ($annot_update_id, $status_link_id, $cdna_acc, $status_id, $is_chromo, $comment, $curated_exception, $curator_comment, $status_descr) = @$result;
        my $ref = {status_link_id => $status_link_id,
                   cdna_acc => $cdna_acc,
                   status_id => $status_id,
                   is_chromo => $is_chromo,
                   comment => $comment,
                   curated_exception => $curated_exception,
                   curator_comment => $curator_comment,
                   status_descr => $status_descr,
                   update_id=>$annot_update_id};
        push (@x, $ref);
    }
    return (@x);
}


sub get_annot_update_info {
    my $dbproc = shift;
    my $update_id = shift;
    
    my $struct;
    my $query = "select * from annotation_updates where update_id = ?";
    my $result = &first_result_sql($dbproc, $query, $update_id);
    if ($result) {
        my ($update_id, $gene_id, $model_id, $alt_splice_flag, $before_gene_obj, $after_gene_obj, $compare_id, $is_valid, $have_before, $have_after) = @$result;
        $struct = {update_id=>$update_id,
                   gene_id=>$gene_id,
                   model_id=>$model_id,
                   alt_splice_flag=>$alt_splice_flag,
                   before_gene_obj=>$before_gene_obj,
                   after_gene_obj=>$after_gene_obj,
                   compare_id=>$compare_id,
                   is_valid=>$is_valid,
                   have_before=>$have_before,
                   have_after=>$have_after};
    }
    return ($struct);
}


####
sub get_fli_status {
    my $dbproc = shift;
    my $acc = shift;
    my $query = "select is_fli from cluster_link where cdna_acc = ?";
    my $result_aref = &first_result_sql($dbproc, $query, $acc);
    return ($result_aref->[0]);
}

####
sub get_annot_db {
    my $dbproc = shift;
    my $query = "select annot_db, annot_server from project";
    my $result = &first_result_sql($dbproc, $query);
    my $db = $result->[0];
    my $server = $result->[1];
    return ($db, $server);
}

####
sub get_max_compare_id {
    my $dbproc = shift;
    my $query = "select max(compare_id) from annotation_compare";
    my $result = &first_result_sql($dbproc, $query);
    return ($result->[0]);
}


sub get_latest_annot_version {
	my ($dbproc) = @_;
	
	## get latest annotation version:
	my $query = "select max(version_id) from annotation_admin";
	my $annot_version = &Mysql_connect::very_first_result_sql($dbproc, $query);
	unless ($annot_version) {
	    confess "Sorry, no version of the annotation to compare to yet.\n";
	}

	return($annot_version);
}




####
sub get_cluster_id_via_subcluster_id {
    my ($dbproc, $subcluster_id) = @_;
    my $query = "select cluster_id from subclusters where subcluster_id = ?";
    my $result = &first_result_sql($dbproc, $query, $subcluster_id);
    return ($result->[0]);
}

####
sub get_updated_gene_obj {
    my ($dbproc, $model_feat_name, $compare_id) = @_;
    my $query = "select after_gene_obj from annotation_updates where model_id = ? and is_valid = 1 and alt_splice_flag = 0 and have_after = 1 and compare_id = ?";
    if (wantarray()) {
        my @blobs;
        my @results = &do_sql_2D($dbproc, $query, $model_feat_name, $compare_id);
        foreach my $result (@results) {
            my $blob = $result->[0];
            if ($blob) {
                push (@blobs, $blob);
            }
        }
        return (@blobs);
    }
    
    my $blob = &very_first_result_sql($dbproc, $query, $model_feat_name, $compare_id);
    return ($blob);
}

sub get_after_gene_via_update_id {
    my ($dbproc, $update_id) = @_;
    my $gene_obj;
    my $query = "select after_gene_obj from annotation_updates where have_after = 1 and update_id = $update_id";
    my $result = &very_first_result_sql($dbproc, $query);
    if ($result) {
        $gene_obj = thaw($result);
    }
    return ($gene_obj);
}

sub get_before_gene {
    my ($dbproc, $model_feat_name, $compare_id) = @_;
    my $query = "select before_gene_obj from annotation_updates where have_before = 1 and compare_id = $compare_id and alt_splice_flag = 0 and model_id = \"$model_feat_name\"";
    my $gene_obj;
    my $blob = &very_first_result_sql($dbproc, $query);
    if ($blob) {
        $gene_obj = thaw($blob);
    }
    return ($gene_obj);
}


####
sub get_seq_from_fasta {
    my ($acc, $fasta_db) = @_;
    unless (-s "$fasta_db.cidx") {
        my $cdbfasta = $ENV{CDBFASTA} || "cdbfasta";
        my $ret = system ("$cdbfasta $fasta_db");
        if ($ret) {
            die "ERROR, couldn't index $fasta_db using cdbfasta.\n";
        }
    }
    
    my $cdbyank = $ENV{CDBYANK} || 'cdbyank';
    
    if ($acc =~ /\.path\d+$/) {
        ## multiple alignments allowed for a single sequence, discriminated by their .path no, GMAP-style
        $acc =~ s/\.path\d+$//;
    }
    

    my $cmd = "$cdbyank -a \'$acc\' $fasta_db.cidx";
    print "CMD: $cmd\n" if $SEE;
    
    my $seq = `$cmd`;
    unless ($seq) { die "FATAL: couldn't retrieve seq for $acc from $fasta_db\n";}
    my @x = split (/\n/, $seq);
    shift @x;
    $seq = join ("", @x);
    $seq =~ s/\s//g;
    return (uc $seq);
}


####
sub get_gene_objs_via_gene_id {
    my ($dbproc, $gene_id, $annot_version) = @_;
    my @gene_objs;
    my $query = "select gene_obj from annotation_store where gene_id = ? and annotation_version = ?";
    my @results = &do_sql_2D($dbproc, $query, $gene_id, $annot_version);
    foreach my $result (@results) {
        my ($gene_obj) = thaw($result->[0]);
        push (@gene_objs, $gene_obj);
    }
    return (@gene_objs);
}

####
sub get_gene_obj_via_model_id {
    my ($dbproc, $model_id, $annot_version) = @_;
    my $query = "select gene_obj from annotation_store where model_id = ? and annotation_version = ?";
    my $gene_obj = thaw(&very_first_result_sql($dbproc, $query, $model_id, $annot_version));
    return ($gene_obj);
}

####
sub get_alignment_span {
    my ($dbproc, $align_acc) = @_;
    
    my $query = "select min(a.lend), max(a.rend) from alignment a, align_link al where al.align_id = a.align_id and al.align_acc = ? and al.validate = 1";
    my $result = &first_result_sql($dbproc, $query, $align_acc);
    return (@$result);
}

####
sub get_asmbl_gene_obj {
    my ($dbproc, $cdna_acc, $allow_5prime_partial, $allow_3prime_partial) = @_;
    
    my $query = "select gene_obj from asmbl_gene_objs where cdna_acc = ? and allow_5prime_partial = ? and allow_3prime_partial = ?";
    my $blob = &very_first_result_sql($dbproc, $query, $cdna_acc, $allow_5prime_partial, $allow_3prime_partial);
    
    my $gene_obj = thaw ($blob);

    unless (ref $gene_obj) {
        confess "Error, no gene_obj returned based on query: $query, $cdna_acc, $allow_5prime_partial, $allow_3prime_partial ";
    }
    
    return ($gene_obj);
}



####
sub load_CDNA_alignment_obj {
    my ($dbproc, $alignment, $prog, $validate) = @_;
    
    my @segments = $alignment->get_alignment_segments();
    my $num_segments = $#segments + 1;
        
    
    ## insert row in cdna_link
    my $query = "insert align_link (align_acc, cdna_info_id, prog, validate, aligned_orient, spliced_orient, num_segments) values (?,?,?,?,?,?,?)";
    &RunMod($dbproc, 
            $query, 
            $alignment->get_acc(),
            $alignment->get_cdna_id(),
            $prog, 
            $validate, 
            $alignment->get_aligned_orientation(), 
            $alignment->get_spliced_orientation(), 
            $num_segments,
            );
    
    my $align_id = &Mysql_connect::get_last_insert_id($dbproc);
        
    foreach my $segment (@segments) {
        my ($lend, $rend) = sort {$a<=>$b} $segment->get_coords();
        my ($mlend, $mrend) = sort {$a<=>$b} $segment->get_mcoords();
        my $query = "insert alignment (align_id, lend, rend, mlend, mrend, orient, per_id) values (?,?,?,?,?,?,?)";
        &RunMod($dbproc, $query, $align_id, $lend, $rend, $mlend, $mrend, $alignment->get_aligned_orientation(), $segment->get_per_id());
        
    }
    
    return ($align_id);
}


1; #EOM.

