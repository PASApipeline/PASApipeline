#!/usr/bin/env perl

use Pasa_init;
use Pasa_conf;
use DBI;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
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
our $IMAGE_X_SIZE = $cgi->param('IMAGE_X_SIZE') || 750;
my $db = $cgi->param('db') or die "Need db";
my $update_id = $cgi->param('update_id') or die "Need update_id";



my $DEBUG = $cgi->param('DEBUG');
if ($DEBUG) {
    print $cgi->header();
    print "<pre>\n";
} else {
    print $cgi->header(-type=>'image/gif');
}


 main: 
{
    my $usage = "usage: $0 subcluster_id";
    
    my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
    my $mysql_ro_user = &Pasa_conf::getParam("MYSQL_RO_USER");
    my $mysql_ro_password = &Pasa_conf::getParam("MYSQL_RO_PASSWORD");

    my ($dbproc) = &connect_to_db($mysql_server,$db,$mysql_ro_user,$mysql_ro_password);
    
    ## Get the assemblies corresponding to the update
    my $query = "select distinct cdna_acc from status_link where annot_update_id = $update_id";
    my @cdnas;
    my @results = &do_sql_2D($dbproc, $query);
    foreach my $result (@results) {
	my $acc = $result->[0];
	my $align_id = &Ath1_cdnas::get_align_id_via_align_acc($dbproc, $acc);
        my $alignment = &Ath1_cdnas::create_alignment_obj($dbproc, $align_id);
        push (@cdnas, $alignment);
    }
    @cdnas = sort {$a->{length}<=>$b->{length}} @cdnas;
    
    ## Get the gene_obj
    my $query = "select compare_id, model_id, after_gene_obj from annotation_updates where update_id = $update_id";
    my $result = &first_result_sql($dbproc, $query);
    if ($result) {
	my ($compare_id, $model_id, $after_blob) = @$result;
	if ($after_blob) {
	    my $after_gene_obj = thaw($after_blob);
	    $after_gene_obj->{Model_feat_name} = "After Update";
	    unshift (@cdnas, $after_gene_obj);
	}

	if ($model_id) {
	    my @model_id_list = split (/\s+/, $model_id);
	    foreach my $model_id (@model_id_list) {
		my $query = "select ast.gene_obj from annotation_store ast, annotation_compare ac where ac.compare_id = $compare_id and ac.annotation_version = ast.annotation_version and ast.model_id = ?";
		my $before_blob = &very_first_result_sql($dbproc, $query, $model_id);
		if ($before_blob) {
		    my $before_gene_obj = thaw($before_blob);
		    $before_gene_obj->{Model_feat_name} .= " Before Update";
		    unshift (@cdnas, $before_gene_obj);
		}
	    }
	}
    


    }    
    
    my $gene_cdna_imager = new Gene_cdna_image;
    
    my $image = $gene_cdna_imager->create_image(@cdnas);
    print $image->png() unless $DEBUG;
}





