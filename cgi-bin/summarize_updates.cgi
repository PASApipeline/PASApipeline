#!/usr/bin/env perl

use Pasa_init;
use Pasa_conf;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Mysql_connect;
use strict;
use DBI;
use Data::Dumper;
use CDNA::CDNA_alignment;
use Gene_obj;
use Ath1_cdnas;
use Getopt::Std;
use Storable qw (thaw);

use vars qw ($opt_h $opt_D $opt_p $opt_d $DEBUG $opt_S $opt_M $opt_X $opt_v);

&getopts ('hD:dp:S:M:Xv');


$|=1;

open (STDERR, ">&STDOUT");

my $cgi = new CGI;
print $cgi->header('text/html');

my $db = $cgi->param('db') or die "Need db";
my $compare_id = $cgi->param('compare_id') or die "need compare_id";
our $SEE = ($cgi->param('SEE') || $cgi->param('DEBUG'));

my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $mysql_ro_user = &Pasa_conf::getParam("MYSQL_RO_USER");
my $mysql_ro_password = &Pasa_conf::getParam("MYSQL_RO_PASSWORD");
my ($dbproc) = &connect_to_db($mysql_server,$db,$mysql_ro_user,$mysql_ro_password);

print $cgi->start_html(-title=>"cdna assembly report");

print "<pre>\n";

my $query = "select distinct a.update_id, a.gene_id, a.model_id, a.alt_splice_flag, a.is_novel_flag "
    . " from annotation_updates a, status s, status_link sl "
    . " where a.compare_id = $compare_id "
    . " and a.is_valid = 1 and a.have_after = 1 "
    . " and sl.compare_id = $compare_id "
    . " and s.status_id = sl.status_id and s.requires_update = 1 "
    . " and sl.annot_update_id = a.update_id";

my @results = &Mysql_connect::do_sql_2D($dbproc, $query);
foreach my $result (@results) {
    my ($update_id, $gene_id, $model_id, $alt_splice_flag, $is_novel_flag) = @$result;
    #print "\n//$update_id, $gene_id, $model_id, altsplice: $alt_splice_flag, novel: $is_novel_flag\n";
    my $query = "select after_gene_obj from annotation_updates where update_id = $update_id";
    my $after_gene_obj = &Mysql_connect::very_first_result_sql($dbproc, $query);
    my $asmbl_id = &get_asmbl_id_via_update_id($update_id);
    #print "asmbl_id: $asmbl_id\n";
    my $new_gene_obj = thaw($after_gene_obj);
    my $comment;
    if ($alt_splice_flag) {
	print "<hr><h2>Add Alternatively Spliced Isoform to Gene_id: $gene_id</h2>\n";
    } elsif ($is_novel_flag) {
	print "<hr><h2>Add a New Gene</h2>\n";
    } elsif ($model_id) { #update
	print "<hr><h2>Update Existing Gene (Gene_id: $gene_id, Model_id: $model_id)</h2>\n";
    } else {
	die "Death: Not sure what to do here.  No actions specified.\n";
    }
    print "Genomic sequence ID: $asmbl_id\n";
    print $new_gene_obj->toString(-noSeqs=>1);
    
}
print "<hr><h2>Updates complete for this round.</h2>\n";


####
sub get_asmbl_id_via_update_id {
    my ($update_id) = @_;

    my $query = "select cdna_acc from status_link where annot_update_id = $update_id";
    my $cdna_acc = &Mysql_connect::very_first_result_sql($dbproc, $query);
    
    my $query = "select c.annotdb_asmbl_id from clusters c, align_link al where al.align_acc = ? and al.cluster_id = c.cluster_id";
    my $asmbl_id = &Mysql_connect::very_first_result_sql($dbproc, $query, $cdna_acc);
    return ($asmbl_id);
}

