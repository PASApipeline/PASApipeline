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


my $status_id = $cgi->param('status_id') or die "Gimme a status id";
my $db = $cgi->param('db') or die "Gimme a db.";
my $compare_id = $cgi->param('compare_id') or die "Need compare_id";
my $annot_db = $cgi->param('annot_db');

my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $mysql_ro_user = &Pasa_conf::getParam("MYSQL_RO_USER");
my $mysql_ro_password = &Pasa_conf::getParam("MYSQL_RO_PASSWORD");
my ($dbproc) = &connect_to_db($mysql_server,$db,$mysql_ro_user,$mysql_ro_password);

print $cgi->start_html(-title=>"assembly listing");
print "<center>\n";

my $query = "select status_descr from status where status_id = ?";
my $result_aref = &first_result_sql($dbproc, $query, $status_id);
my $status_descr = $result_aref->[0];

print "<h1>Status($status_id): $status_descr</h1>\n";

my $query = "select subl.subcluster_id, subl.cdna_acc from subcluster_link subl, status_link sl where subl.cdna_acc = sl.cdna_acc and sl.status_id = $status_id and sl.compare_id = $compare_id";

my %offending_asmbl; #store offending asmbl_info.
my @results = &do_sql_2D($dbproc, $query);
my %data; #keyed by subcluster_id
my %subclusters;
foreach my $result_aref (@results) {
    my ($subcluster_id, $cdna_acc) = @$result_aref;
    $offending_asmbl{$cdna_acc} = 1;
    $subclusters{$subcluster_id} = 1;
}

foreach my $subcluster_id (keys %subclusters) {
    my $query = "select subl.cdna_acc from subcluster_link subl where subl.subcluster_id = ?";
    my @results = &do_sql_2D($dbproc, $query, $subcluster_id);
    my @cdnas;
    foreach my $result_aref (@results) {
	my $cdna = $result_aref->[0];
	push (@cdnas, $cdna);
    }
    $data{$subcluster_id} = [@cdnas];
}

my $i = 0;
print "<table border=3 align=center><tr><th>#</th><th>subcluster id</th><th>assembly acc</th><th>Chromosome Annotations</th></tr>\n";
foreach my $subcluster_id ( sort {$#{$data{$a}}<=>$#{$data{$b}} || $a<=>$b} keys %data) {
    my @cdnas = @{$data{$subcluster_id}};
    $i++;
    print "<tr><td><b>$i</b></td><td><a href = \"cdnaasmbl_report.cgi?subcluster_id=$subcluster_id&db=$db&compare_id=$compare_id\">$subcluster_id</a></td><td>";
    print "<table>\n";
    foreach my $cdna (@cdnas) {
	my $bgcolor = "#FFFFFF";
	if ($offending_asmbl{$cdna}) {
	    $bgcolor = "#FF0000";
	}
	print "<tr><td bgcolor = \"$bgcolor\">$cdna</td></tr>\n";
    }
    print "</table>\n</td><td width=\"50%\">";
    
    my %genes;
    foreach my $cdna (@cdnas) { 
	my $query = "select gene_id from annotation_link where cdna_acc = ? and compare_id = $compare_id";
	my @results = &do_sql_2D($dbproc, $query, $cdna);
	
	foreach my $result (@results) {
	    my $gene_id = $result->[0];
	    $genes{$gene_id} = 1;
	}
    }
    my @genes = keys %genes;
    
    print "@genes</td></tr>\n";
}
print "</table>\n";

print $cgi->end_html();

$dbproc->disconnect;












