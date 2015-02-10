#!/usr/bin/env perl

use Pasa_init;
use Pasa_conf;
use Mysql_connect;
use strict;
use DBI;
use Data::Dumper;
use Getopt::Std;
use Ath1_cdnas;
use CDNA::CDNA_alignment;
use CGI;
use CGI::Carp qw(fatalsToBrowser);


open (STDERR, "&>STDOUT");


my $cgi = new CGI();

my $db = $cgi->param('db') or die "Need db param\n";

print $cgi->header('text/plain'); #now prompts to save.

my ($MYSQLserver, $MYSQLdb) = ("localhost", $db);
my ($user, $password) = ("access", "access");

my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $mysql_ro_user = &Pasa_conf::getParam("MYSQL_RO_USER");
my $mysql_ro_password = &Pasa_conf::getParam("MYSQL_RO_PASSWORD");
my ($dbproc) = &connect_to_db($mysql_server,$db,$mysql_ro_user,$mysql_ro_password);

my $query = "select accession, sequence from cdna_sequence";
my @results = &do_sql_2D($dbproc, $query);
foreach my $result (@results) {
    my ($accession, $sequence) = @$result;
    $sequence =~ s/(\w{60})/$1\n/g;
    print ">$accession.$db\n$sequence\n";
}


$dbproc->disconnect;



	
