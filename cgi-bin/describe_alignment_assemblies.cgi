#!/usr/bin/env perl

use Pasa_init;
use Pasa_conf;
use Mysql_connect;
use strict;
use DBI;
use Data::Dumper;
use Gene_obj;
use Getopt::Std;
use Storable qw (freeze thaw);
use CGI;
use CGI::Carp qw(fatalsToBrowser);


my $cgi = new CGI;
print $cgi->header('text/plain');
my $db = $cgi->param('db');


$|=1;
our $SEE = 0;

open (STDERR, "&>STDOUT");


my ($MYSQLdb, $MYSQLserver) = ($db,"localhost");
my $INCLUDE_ASSEMBLY_ALIGNMENT = 1;
#process passwords
my ($user, $password) = ("access", "access");

my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $mysql_ro_user = &Pasa_conf::getParam("MYSQL_RO_USER");
my $mysql_ro_password = &Pasa_conf::getParam("MYSQL_RO_PASSWORD");
my ($dbproc) = &connect_to_db($mysql_server,$db,$mysql_ro_user,$mysql_ro_password);


my %assembly_to_asmbl_id;

my $query = "select c.annotdb_asmbl_id, al.align_acc from clusters c, align_link al, cdna_info ci "
    . " where c.cluster_id = al.cluster_id and al.cdna_info_id = ci.id and ci.is_assembly = 1";
my @results = &Mysql_connect::do_sql_2D($dbproc, $query);
foreach my $result (@results) {
    my ($asmbl_id, $cdna_acc) = @$result;
    $assembly_to_asmbl_id{$cdna_acc} = $asmbl_id;
}

## Get listing of cdna_accs for each asmbl_id:
my %asmbl_to_accs;
my $query = "select asmbl_acc, cdna_acc from asmbl_link";
my @results = &Mysql_connect::do_sql_2D($dbproc, $query);
foreach my $result (@results) {
    my ($asmbl_acc, $cdna_acc) = @$result;
    my $listref = $asmbl_to_accs{$asmbl_acc};
    unless (ref $listref) {
	$listref = $asmbl_to_accs{$asmbl_acc} = [];
    }
    push (@$listref, $cdna_acc);
}

foreach my $asmbl (sort keys %asmbl_to_accs) {
    my $listref = $asmbl_to_accs{$asmbl};
    my $cdna_accs = join (",", @$listref);
    my $asmbl_id = $assembly_to_asmbl_id{$asmbl};
    my $assembly_alignment = "";
    if ($INCLUDE_ASSEMBLY_ALIGNMENT) {
	my $query = "select alignment from align_link where align_acc = \"$asmbl\"";
	$assembly_alignment = &Mysql_connect::very_first_result_sql($dbproc, $query);
    }
    if ($assembly_alignment) {
	$assembly_alignment = "alignment_assembly: " . $assembly_alignment;
    }
    print "alignment_assembly:$asmbl\tcdna_accs:$cdna_accs\tgenomic_acc:${asmbl_id}\t${assembly_alignment}\n";
}

$dbproc->disconnect;

exit(0);

