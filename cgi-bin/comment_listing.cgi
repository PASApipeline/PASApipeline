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

$|++;


my $cgi = new CGI;
print $cgi->header('text/html');
my $db = $cgi->param('db');
our $SEE = $cgi->param('SEE') || $cgi->param('DEBUG');
my $compare_id = $cgi->param('compare_id');
my $status_id = $cgi->param('status_id');


unless ($db && defined($status_id) && defined($compare_id)) {
    die "Must set the db, status_id, compare_id parameters.\n";
}

my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $mysql_ro_user = &Pasa_conf::getParam("MYSQL_RO_USER");
my $mysql_ro_password = &Pasa_conf::getParam("MYSQL_RO_PASSWORD");
my ($dbproc) = &connect_to_db($mysql_server,$db,$mysql_ro_user,$mysql_ro_password);
print $cgi->start_html();

my $query = "select status_id, status_descr from status";
my @results = &do_sql_2D($dbproc, $query);

my %status;
foreach my $result (@results) {
    my ($status_id, $status_descr) = @$result;
    $status{$status_id} = $status_descr;
}

my %subclusters;
my $query = "select cdna_acc, subcluster_id from subcluster_link";
my @results = &do_sql_2D($dbproc,$query);
foreach my $result (@results) {
    my ($cdna_acc, $subcluster_id) = @$result;
    $subclusters{$cdna_acc} = $subcluster_id;
}

print "<h1>Comment listing for status_id: $status_id ($status{$status_id})</h1>";

my $query = "select cdna_acc, comment from status_link where status_id = $status_id and compare_id = $compare_id order by comment";
my @results = &do_sql_2D($dbproc, $query);
my $num = 0;
print "<table border=2><tr><th>#</th><th>assembly</th><th>comment</th></tr>\n";

foreach my $result (@results) {
    $num++;
    my ($cdna_acc, $comment) = @$result;
    my $subcluster_id = $subclusters{$cdna_acc};
    print "<tr><td>$num</td><td><a href=\"cdnaasmbl_report.cgi?subcluster_id=$subcluster_id&db=$db&compare_id=$compare_id\">$cdna_acc</a></td><td>$comment</td></tr>\n";
}
print "</table>\n";

print $cgi->end_html();


exit(0);
