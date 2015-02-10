#!/usr/bin/env perl

use Pasa_init;
use Pasa_conf;
use DBI;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use lib ($ENV{EGC_SCRIPTS}, $ENV{MYSQL_LIB}, $ENV{EUK_MODULES});
use Data::Dumper;
use Mysql_connect;
use Ath1_cdnas;
use strict;
use CGI::Pretty ":standard";
use Pasa_CGI;

$|++;

my $cgi = new CGI;
print $cgi->header('text/html');
my $db = $cgi->param('db');
our $SEE = $cgi->param('SEE') || $cgi->param('DEBUG');


my $css_common_text = &Pasa_CGI::get_common_CSS();
print $cgi->start_html(-title=>"PASA assembly clusters for $db",
                       -head => style( { type => "text/css" }, $css_common_text ),
                       );



my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $mysql_ro_user = &Pasa_conf::getParam("MYSQL_RO_USER");
my $mysql_ro_password = &Pasa_conf::getParam("MYSQL_RO_PASSWORD");
my ($dbproc) = &connect_to_db($mysql_server,$db,$mysql_ro_user,$mysql_ro_password);





## Get max compare_id
my $query = "select compare_id from annotation_compare order by compare_id desc";
my $compare_id = &very_first_result_sql($dbproc, $query) || 1;

print "<center><h1>Subcluster summary for $db\n";



my $query = "select subcluster_id, count(*) from subcluster_link group by subcluster_id order by 2 desc";
my @results = &do_sql_2D($dbproc, $query);

my $outTable = "<table border=1 ><tr><th>Subcluster_ID</th><th># pasa assemblies</th></tr>\n";
my %subclusterSizeCounter;

my $counter = 0;
foreach my $result (@results) {
    my ($subcluster_id, $count) = @$result;
    $counter++;
    $subclusterSizeCounter{$count}++;
    $outTable .= "<tr><td align=center ><a href = \"cdnaasmbl_report.cgi?subcluster_id=$subcluster_id&db=$db&compare_id=$compare_id\" target=_blank >$subcluster_id</a></td><td align=center >$count</td></tr>\n";
}
$outTable .= "</table>\n";


print "<h2>Subcluster size distribution</h2>\n";
print "<table border=1 ><tr><th>subcluster size</th><th># subclusters</th></tr>\n";
foreach my $size (reverse sort {$a<=>$b} keys %subclusterSizeCounter) {
    print "<tr><td align=center >$size</td><td align=center >$subclusterSizeCounter{$size}</td></tr>\n";
}
print "</table>\n";
print "<h2>Individual Subclusters</h2>\n";
print $outTable;

print $cgi->end_html();

exit(0);


