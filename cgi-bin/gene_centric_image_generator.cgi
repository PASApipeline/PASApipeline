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

my $TU_feat_name = $cgi->param('TU');
our $IMAGE_X_SIZE = $cgi->param('IMAGE_X_SIZE') || 750;
my $db = $cgi->param('db');
my $compare_id = $cgi->param('compare_id') or die "Need compare_id";
unless ($compare_id && $TU_feat_name && $db) {
    die "Need db, compare_id and TU_feat_name";
}


my $DEBUG = $cgi->param('DEBUG');
if ($DEBUG) {
    print $cgi->header();
    print "<pre>\n";
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
    
    # get annot version
    my $query = "select annotation_version from annotation_compare where compare_id = ?";
    my $annot_version = &very_first_result_sql($dbproc, $query, $compare_id);
	
	
    # show the gene model in annot_db, the updated forms, and all corresponding assemblies.
    my $query = "select sl.cdna_acc, sl.annot_update_id "
		. "from status_link sl, annotation_link al, status s "
		. "where al.gene_id = ? and al.cdna_acc = sl.cdna_acc and al.compare_id = $compare_id and sl.compare_id = $compare_id"
		. " and s.status_id = sl.status_id and s.fails_incorporation = 0"
		
		;

    my @results = &do_sql_2D($dbproc, $query, $TU_feat_name);
    my %update_ids;
    my %cdna_accs;
    foreach my $result (@results) {
		my ($cdna_acc, $update_id) = @$result;
		$cdna_accs{$cdna_acc} = 1;
		if ($update_id) {
			$update_ids{$update_id} = 1;
		}
    }
    
    my @geneobjs;
    ## Get current incarnation of gene:
    my @genes = &Ath1_cdnas::get_gene_objs_via_gene_id($dbproc, $TU_feat_name, $annot_version);
    my @model_feats;
    foreach my $gene (@genes) {
		$gene->{com_name} = "current structure";
		push (@model_feats, $gene->{Model_feat_name});
    }
    push (@geneobjs, @genes);
    
    ## Get updates:
    foreach my $model_feat (@model_feats) {
		my (@blobs) = &Ath1_cdnas::get_updated_gene_obj($dbproc, $model_feat, $compare_id);
		if (@blobs) {
			foreach my $blob (@blobs) {
				my $gene = thaw($blob);
				$gene->{com_name} = "tentative update";
				push (@geneobjs, $gene);
			}
		}
    }
    
    ## Get the alt splice isoforms:
    my $query = "select after_gene_obj from annotation_updates where compare_id = $compare_id and gene_id = \"$TU_feat_name\" and alt_splice_flag = 1 and is_valid = 1 and have_after = 1";
    my @results = &Mysql_connect::do_sql_2D($dbproc, $query);
    foreach my $result (@results) {
		my $blob = $result->[0];
		if ($blob) {
			my $gene = thaw($blob);
			$gene->{com_name} = "tentative alt-splicing isoform.";
			push (@geneobjs, $gene);
		}
    }
    
	
    ## Get failed updates:
    ## See if gene structures are available for each assembly
    foreach my $update_id (keys %update_ids) {
		my $query = "select after_gene_obj from annotation_updates where update_id = ? and is_valid = 0";
		my $result = &Mysql_connect::very_first_result_sql($dbproc, $query, $update_id);
		if ($result) {
			my $gene_obj = thaw($result);
			$gene_obj->{com_name} = "cdna failed update"; 
			print $gene_obj->toString() if $DEBUG;
			push (@geneobjs, $gene_obj);
		}
    }
    
    ## Get the assembly cdna objects:
    foreach my $cdna_acc (keys %cdna_accs) {
		my $align_id = &Ath1_cdnas::get_align_id_via_align_acc($dbproc, $cdna_acc);
		my $alignment = &Ath1_cdnas::create_alignment_obj($dbproc, $align_id);
		
		push (@geneobjs, $alignment);
    }
	
    my $gene_cdna_imager = new Gene_cdna_image;
    
    my $image = $gene_cdna_imager->create_image(@geneobjs);
    print $image->png() unless $DEBUG;
    
    $dbproc->disconnect;
    
}

exit(0);





