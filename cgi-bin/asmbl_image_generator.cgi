#!/usr/bin/env perl

use FindBin;
use lib ($FindBin::Bin);
use Pasa_init;
use Pasa_conf;
use DBI;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Sequence_coords_image;
use Data::Dumper;
use Gene_obj;
use CDNA::CDNA_alignment;
use Mysql_connect;
use Ath1_cdnas;
use strict;
use Storable qw (thaw);
use Gene_cdna_image;

$| = 1;

my $cgi = new CGI;


my $subcluster_id = $cgi->param('subcluster_id');
my $assembly_acc = $cgi->param('assembly_acc');
my $IMAGE_X_SIZE = $cgi->param('IMAGE_X_SIZE') || 750;
my $DRAW_PANEL_SCALER = $cgi->param('DRAW_PANEL_SCALER') || 0.7;

my $db = $cgi->param('db');
my $compare_id = $cgi->param('compare_id') or die "Need compare_id";
unless ($subcluster_id or $assembly_acc) {
    die "Need subcluster_id or assembly_acc";
}


my $DEBUG = $cgi->param('DEBUG');
if ($DEBUG) {
    our $SEE = 1;
    print $cgi->header();
    print "<pre>\n";
    print "INC: @INC\nPerllib: $ENV{PERLLIB}\n";
} else {
    print $cgi->header(-type=>'image/gif');
}

my $consensus = 10;
my $consensus_color = "blue";

 main: 
{
    my $usage = "usage: $0 subcluster_id";
    
    my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
    my $mysql_ro_user = &Pasa_conf::getParam("MYSQL_RO_USER");
    my $mysql_ro_password = &Pasa_conf::getParam("MYSQL_RO_PASSWORD");
    my ($dbproc) = &connect_to_db($mysql_server,$db,$mysql_ro_user,$mysql_ro_password);
    my (%TUs, %assemblies, %cdnas);
    unless ($assembly_acc xor $subcluster_id) { die "Use assembly_acc or subcluster_id as parameters.\n";}
    if ($subcluster_id) {
        ## From the subcluster_id, get the genes linked, assembly accs, and the sub_accs:
        my $query = "select cdna_acc from subcluster_link where subcluster_id = ?";
        my @results = &Mysql_connect::do_sql_2D ($dbproc, $query, $subcluster_id);
        
        foreach my $result (@results) {
            my $assembly_acc = $result->[0];
            $assemblies{$assembly_acc} = 1;
        }
    } else {
        $assemblies{$assembly_acc} = 1;
    }
    
    ## Get the annotation version:
    my $query = "select annotation_version from annotation_compare where compare_id = ?";
    my $annot_version = &Mysql_connect::very_first_result_sql($dbproc, $query, $compare_id);
    
    ## Get the genes linked to and cDNAs built into each assembly.
    foreach my $assembly_acc (keys %assemblies) {
        
        ## Get the genes linked:
        my $query = "select gene_id from annotation_link where cdna_acc = ? and compare_id = ?";
        my @results = &Mysql_connect::do_sql_2D ($dbproc, $query, $assembly_acc, $compare_id);
        foreach my $result (@results) {
            my $tu = $result->[0];
            $TUs{$tu} = 1;
        }
        
        ## Get the cdnas comprising the assembly
        my $query = "select cdna_acc from asmbl_link where asmbl_acc = ?";
        my @results = &Mysql_connect::do_sql_2D ($dbproc, $query, $assembly_acc);
        foreach my $result (@results) {
            my $cdna_acc = $result->[0];
            $cdnas{$cdna_acc} = 1;
        }
    }
    
    ## convert to structs of data:
    my @objs;
    foreach my $tu (keys %TUs) {
        my $query = "select model_id from annotation_store where gene_id = ? and annotation_version = ?";
        my @results = &do_sql_2D($dbproc, $query, $tu, $annot_version);
        my @model_feat_names;
        foreach my $result (@results) {
            push (@model_feat_names, $result->[0]);
        }
        
        foreach my $model_feat_name (@model_feat_names) {
            my $gene_obj = Ath1_cdnas::get_before_gene($dbproc, $model_feat_name, $compare_id);
            if (ref $gene_obj) {
                $gene_obj->{com_name} = "[before]: " . $gene_obj->{com_name};
                push (@objs, $gene_obj);
            }
            $gene_obj = Ath1_cdnas::get_gene_obj_via_model_id($dbproc, $model_feat_name, $annot_version);
            $gene_obj->{com_name} = "[current(v$annot_version)]: " . $gene_obj->{com_name};
            push (@objs, $gene_obj);
        }
    }
    
    ## See if gene structures are available for each assembly
    foreach my $assembly_acc (keys %assemblies) {
        my $query = "select annot_update_id from status_link where compare_id = ? and cdna_acc = ?";
        my $result = &Mysql_connect::first_result_sql($dbproc, $query, $compare_id, $assembly_acc);
        if ($result) {
            my $annot_update_id = $result->[0];
            print "Asmbl_acc: $assembly_acc, Annot_update_id: $annot_update_id\n" if $DEBUG;
            my $query = "select after_gene_obj from annotation_updates where update_id = ?";
            my $result = &Mysql_connect::first_result_sql($dbproc, $query, $annot_update_id);
            if ((ref $result) && $result->[0]) {
                my $blob = $result->[0];
                my $gene_obj = thaw($blob);
                $gene_obj->{Model_feat_name} = "$assembly_acc-including gene model"; 
                print $gene_obj->toString() if $DEBUG;
                push (@objs, $gene_obj);
            }
        }
    }
    
    my @cdna_data;
    foreach my $acc (keys %assemblies, keys %cdnas) {
        #my $align_id = &Ath1_cdnas::get_validating_align_id_via_acc($dbproc, $acc);
        #my $alignment = &Ath1_cdnas::create_alignment_obj($dbproc, $align_id);
       
        my $alignment = &Ath1_cdnas::get_alignment_obj_via_align_acc($dbproc, $acc);

        my $orientation = $alignment->get_orientation();
        if ($orientation !~ /[\+\-]/) {
            die "Error.";
        }

        my $query = "select genomicCoord from transcriptPolyA t, align_link al where al.align_acc = ? and al.align_id = t.align_id";
        my $genomicCoord = &very_first_result_sql($dbproc, $query, $acc);
        if ($genomicCoord) {
            $alignment->{polyAcoord} = $genomicCoord; # abusing the alignmentObj here
        }
        
        push (@cdna_data, $alignment);
    }
    @cdna_data = reverse sort {$a->{length}<=>$b->{length}} @cdna_data;
    push (@objs, @cdna_data);
    
    my $gene_cdna_imager = new Gene_cdna_image;
    $gene_cdna_imager->{Sequence_feature_drawer_obj}->{IMAGE_X_SIZE} = $IMAGE_X_SIZE;
    $gene_cdna_imager->{Sequence_feature_drawer_obj}->{DRAW_PANEL_SCALER} = $DRAW_PANEL_SCALER;
    
    
    my $image = $gene_cdna_imager->create_image(@objs);
    print $image->png() unless $DEBUG;
}





