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
my $no_cache = $cgi->param('NOCACHE');
my $compare_id = $cgi->param('compare_id');


unless ($db && $compare_id) {
    die "Must set the db and compare_id parameters.\n";
}

my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $mysql_ro_user = &Pasa_conf::getParam("MYSQL_RO_USER");
my $mysql_ro_password = &Pasa_conf::getParam("MYSQL_RO_PASSWORD");
my ($dbproc) = &connect_to_db($mysql_server,$db,$mysql_ro_user,$mysql_ro_password);

my @models_supported_by_PASA_assemblies;
my @models_supported_by_FL_assemblies;
my @models_supported_by_EST_assemblies_only;


my @models_failed_by_PASA_assemblies;
my @models_failed_by_FL_assemblies;
my @models_failed_by_EST_assemblies_only;


################# Analyze Successes ##############################

my $root_query = "select distinct gene_id, model_id "
    . " from annotation_link al, status_link sl, status s, align_link al2, cdna_info ci "
    . " where al.cdna_acc = sl.cdna_acc "
    . " and sl.status_id = s.status_id "
    . " and sl.compare_id = $compare_id "
    . " and al.compare_id = $compare_id "
    . " and sl.cdna_acc = al2.align_acc "
    . " and al2.cdna_info_id = ci.id "
    . " and ci.is_assembly = 1 ";

my $query = $root_query;
$query .= "and (s.requires_update = 1 or (s.requires_update = 0 and s.fails_incorporation = 0)) ";
my @results = &do_sql_2D($dbproc, $query);
foreach my $result (@results) {
    my $model = join (", ", @$result);
    push (@models_supported_by_PASA_assemblies , $model);
}

my %fl;
## get those considered to be FL-only:
$query .= " and ci.is_fli = 1";
my @results = &do_sql_2D($dbproc, $query);
foreach my $result (@results) {
    my $model = join (", ", @$result);
    $fl{$model} = 1;
    push (@models_supported_by_FL_assemblies, $model);
}

## find those that are not FL:
foreach my $model (@models_supported_by_PASA_assemblies) {
    unless ($fl{$model}) {
        push (@models_supported_by_EST_assemblies_only, $model);
    }
}



################## Analyze Failures #################################
$query = $root_query;
$query .= " and s.fails_incorporation = 1";
my @results = &do_sql_2D($dbproc, $query);
foreach my $result (@results) {
    my $model = join (", ", @$result);
    push (@models_failed_by_PASA_assemblies, $model);
}

## get those only FL
%fl = (); #reinit
$query .= " and ci.is_fli = 1";
my @results = &do_sql_2D($dbproc, $query);
foreach my $result (@results) {
    my $model = join (", ", @$result);
    $fl{$model} = 1;
    push (@models_failed_by_FL_assemblies, $model);
}

## get those failed w/ ESTs only:
foreach my $model (@models_failed_by_PASA_assemblies) {
    unless ($fl{$model}) {
	push (@models_failed_by_EST_assemblies_only, $model);
    }
}


##########################################
#### Summarize DATA ######################
##########################################

print "<p>For successful incorporations, all models (isoforms) of a given gene incorporating alignment assemblies are reported separately.  In the case of failures, no model of a given gene was capable of incorporating the alignment assembly; only one such model is reported for example purposes.  The total number of models needed to incorporate all failed transcripts is unknown.</p>";

print "<table border=1 >\n"
    . "<tr><th colspan=2 >Successful Incorporations</th></tr>\n"
    . "<tr><td>models supported by PASA assemblies</td><td align=right >" . scalar (@models_supported_by_PASA_assemblies) . "</td></tr>\n"
    . "<tr><td>models supported by FL assemblies</td><td align=right >" . scalar (@models_supported_by_FL_assemblies) . "</td></tr>\n"
    . "<tr><td>models supported by EST assemblies only</td><td align=right>" . scalar (@models_supported_by_EST_assemblies_only) . "</td></tr>\n"
    . "<tr><th colspan=2>Failed Incorporations</th></tr>\n"
    . "<tr><td>models failing PASA assembly incorporation</td><td align=right>" . scalar (@models_failed_by_PASA_assemblies) . "</td></tr>\n"
    . "<tr><td>models failing FL assembly incorporation</td><td align=right>" . scalar (@models_failed_by_FL_assemblies) . "</td></tr>\n"
    . "<tr><td>models failing EST assemblies only</td><td align=right>" . scalar (@models_failed_by_EST_assemblies_only) . "</td></tr>\n"
    . "</table>\n";

    

print "<p><pre>";


############### Print out raw data ############################

print "// All Models supported by PASA assemblies:\n";
&dump_data(@models_supported_by_PASA_assemblies);

print "\n\n// Models supported by FL assemblies only:\n";
&dump_data(@models_supported_by_FL_assemblies);

print "\n\n// Models supported by EST assemblies only:\n";
&dump_data(@models_supported_by_EST_assemblies_only);

print "\n\n// Models failed by PASA assemblies:\n";
&dump_data(@models_failed_by_PASA_assemblies);

print "\n\n// Models failed by FL assemblies:\n";
&dump_data(@models_failed_by_FL_assemblies);

print "\n\n// Models failed by EST assemblies only:\n";
&dump_data(@models_failed_by_EST_assemblies_only);

print "</pre>";

print $cgi->end_html();

exit(0);


####
sub dump_data {
    my @models = @_;
    foreach my $model (@models) {
	print "$model\n";
    }
}







