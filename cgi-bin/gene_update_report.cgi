#!/usr/bin/env perl

use Pasa_init;
use Pasa_conf;
use DBI;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Data::Dumper;
use Mysql_connect;
use Ath1_cdnas;
use strict;

my $cgi = new CGI;
print $cgi->header('text/html');

my $step = 10;
my $db = $cgi->param('db') or die "Need db";
my $compare_id = $cgi->param('compare_id') or die "need compare_id";
my $status_id = $cgi->param('status_id') or die "need status_id";
my $count = $cgi->param('count') || 0;
my $stop = $count + $step;
my $annot_db = $cgi->param('annot_db');

my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $mysql_ro_user = &Pasa_conf::getParam("MYSQL_RO_USER");
my $mysql_ro_password = &Pasa_conf::getParam("MYSQL_RO_PASSWORD");
my ($dbproc) = &connect_to_db($mysql_server,$db,$mysql_ro_user,$mysql_ro_password);

print $cgi->start_html(-title=>"gene update report");

my $query = "select status_id, status_descr from status";
my @results = &do_sql_2D($dbproc, $query);
my %status;
foreach my $result (@results) {
    my ($status_id, $status_descr) = @$result;
    $status{$status_id} = $status_descr;
}

print "<h1>Annotation Updates: $status{$status_id} </h1>\n";

my $query = "select distinct au.update_id, au.gene_id, au.model_id from annotation_updates au, status_link sl where sl.compare_id = $compare_id and sl.status_id = $status_id and sl.annot_update_id = au.update_id and au.is_valid = 1 order by au.update_id";
my @results = &do_sql_2D($dbproc, $query);
print "<table border = 2><tr><th>Update ID</th><th>Gene</th><th>Model</th><th>assemblies</th><th>Illustration</th></tr>\n";

my $x = 0;
foreach my $result (@results) {
    $x++;
    if ($x > $count) {
	my ($update_id, $gene_id, $model_id) = @$result;
	## Get the subcluster_id and assembly accessions corresponding to this update:
	my $query = "select distinct sl.subcluster_id, sl.cdna_acc from subcluster_link sl, status_link s where sl.cdna_acc = s.cdna_acc and s.compare_id = $compare_id and s.annot_update_id = $update_id";
	my @results = &do_sql_2D($dbproc, $query);
	my $assembly_text = "<table>";
	foreach my $result (@results) {
	    my ($subcluster_id, $assembly_acc) = @$result;
	    $assembly_text .= "<tr><td><a href=\"cdnaasmbl_report.cgi?db=$db&subcluster_id=$subcluster_id&compare_id=$compare_id\">$assembly_acc</a></td></tr>";
	}
	$assembly_text .= "</table>";
	print "<tr><td>$update_id</td><td>$gene_id<a href=\"gene_centric_image_generator.cgi?db=$db&compare_id=$compare_id&TU=$gene_id\">AllGeneUpd</a></td><td>$model_id</td><td>$assembly_text</td><td><img src=\"annot_update_image.cgi?db=$db&update_id=$update_id\" alt=update_image_$gene_id ></td></tr>\n";
    }
    if ($x > $stop) {
	last;
    }
}

print "</table>\n";
print "<p><a href = \"gene_update_report.cgi?db=$db&count=$x&compare_id=$compare_id&status_id=$status_id\">Next $step results &gt;&gt;</a></p>";

print $cgi->end_html();

$dbproc->disconnect;
exit(0);

