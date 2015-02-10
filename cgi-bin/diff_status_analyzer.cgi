#!/usr/bin/env perl

use strict;
use Pasa_init;
use Pasa_conf;
use DBI;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Data::Dumper;
use Mysql_connect;
use Ath1_cdnas;

$|++;

my $cgi = new CGI;
print $cgi->header('text/html');
my $db = $cgi->param('db');
our $SEE = $cgi->param('SEE') || $cgi->param('DEBUG');

my $compare_id_A = $cgi->param('compare_id_A');
my $compare_id_B = $cgi->param('compare_id_B');

unless ($compare_id_A && $compare_id_B) {
    die "Must specify pair of comparisons via compare_id_A and compare_id_B\n";
}

unless ($db) {
    die "Must set the db parameter.\n";
}

unless ($ENV{WEBSERVER_TMP}) { $ENV{WEBSERVER_TMP} = "/tmp";}

if ($SEE) {
    print "<PRE>\n";
}


my $cache_file = "$ENV{WEBSERVER_TMP}/$db.status.html";
my $time1 = time();
print $cgi->start_html(-title=>"Diff between comparisons: $compare_id_A and $compare_id_B");

my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $mysql_ro_user = &Pasa_conf::getParam("MYSQL_RO_USER");
my $mysql_ro_password = &Pasa_conf::getParam("MYSQL_RO_PASSWORD");
my ($dbproc) = &connect_to_db($mysql_server,$db,$mysql_ro_user,$mysql_ro_password);
    
my %requires_update;
my $query = "select status_id, requires_update from status";
my @results = &do_sql_2D($dbproc, $query);
foreach my $result (@results) {
    my ($status_id, $requires_update) = @$result;
    $requires_update{$status_id} = $requires_update;
}

my $text .= "<hr>\n<h1>Diffs between cDNA Assembly Annotation Comparisons $compare_id_A and $compare_id_B</h1>\n";

my %status_info;

my $query = "select s.status_id, s.requires_update, s.fails_incorporation, s.status_descr, s.rank from status s order by s.rank";
my @results = &do_sql_2D($dbproc, $query);
foreach my $result (@results) {
    my ($status_id, $requires_update, $fails_incorporation, $status_descr, $rank) = @$result;
    $status_info{$status_id} = { requires_update => $requires_update,
				 fails_incorporation => $fails_incorporation,
				 status_descr => $status_descr,
				 rank => $rank };
}


my %Status_A_and_B;

foreach my $compare_id ($compare_id_A, $compare_id_B) {
    my $query = "select cdna_acc, status_id from status_link where compare_id = $compare_id";
    my @results = &do_sql_2D($dbproc, $query);
    foreach my $result (@results) {
	my ($cdna_acc, $status_id) = @$result;
	$Status_A_and_B{$compare_id}->{$cdna_acc} = $status_id;
    }
}

my %cdnas_shifted;
my $A_cdna_set_href = $Status_A_and_B{$compare_id_A};
my $B_cdna_set_href = $Status_A_and_B{$compare_id_B};

foreach my $cdna (keys %$A_cdna_set_href) {
    my $status_id_A = $A_cdna_set_href->{$cdna};
    my $status_id_B = $B_cdna_set_href->{$cdna};
    
    if ($status_id_A != $status_id_B) {
	my $list_aref = $cdnas_shifted{$status_id_A};
	unless (ref $list_aref) {
	    $list_aref = $cdnas_shifted{$status_id_A} = [];
	}
	
	push (@$list_aref, { cdna => $cdna, new_status_id => $status_id_B } );
    }
}


my $diff_text = "";


$text .= "<table border=3 align=center><tr><th>status id</th><th>status description</th><th>Requires Update</th><th>Fails Incorporation</th><th>Number of cDNA assemblies</th><th>Comment Listing</th></tr>\n";
foreach my $status_id (sort {$status_info{$a}->{rank} <=> $status_info{$b}->{rank}} keys %status_info) {
    my $status_struct = $status_info{$status_id};
    
    my ($requires_update, $fails_incorporation, $status_descr) = ($status_struct->{requires_update},
								  $status_struct->{fails_incorporation},
								  $status_struct->{status_descr});
    
    my $update_bgcolor = ($requires_update) ? "bgcolor=\"#00FF00\"" : "";
    my $incorp_bgcolor = ($fails_incorporation) ? "bgcolor=\"#FF0000\"" : "";
    
    my $shifted_entries_aref = $cdnas_shifted{$status_id} || [];
    
    my $count = scalar (@$shifted_entries_aref);
    
    $text .= "<tr><td>$status_id</td><td>$status_descr</td><td $update_bgcolor>&nbsp;</td><td $incorp_bgcolor>&nbsp;</td><td>$count</td><td><a href=\"comment_listing.cgi?db=$db&status_id=$status_id&compare_id=$compare_id_A\">Comments</a></td></tr>\n";
    
    if ($count) {
	&add_diff_text ($status_id, $shifted_entries_aref);
    }
}
$text .= "</table><p>&nbsp;<p>&nbsp;\n";

print $text;
print $diff_text;

print $cgi->end_html();


exit(0);




####
sub add_diff_text {
    my ($orig_status_id, $shifted_entries_aref) = @_;
    
    my %status_id_to_acc_list;
    foreach my $entry (@$shifted_entries_aref) {
	my ($cdna, $new_status_id) = ($entry->{cdna}, $entry->{new_status_id});
	
	my $list_aref = $status_id_to_acc_list{$new_status_id};
	unless (ref $list_aref) {
	    $list_aref = $status_id_to_acc_list{$new_status_id} = [];
	}

	push (@$list_aref, $cdna);
    }
    
    #print "<pre>" . Dumper (\%status_id_to_acc_list);
    
    $diff_text .= "<h2>Changes from Status($orig_status_id) " . $status_info{$orig_status_id}->{status_descr} . "</h2>\n";
    $diff_text .= "<table border=3 align=center><tr><th>new status id</th><th>status description</th><th>Requires Update</th><th>Fails Incorporation</th><th>Accession Listing</th></tr>\n";
    foreach my $status_id (sort {$status_info{$a}->{rank}<=>$status_info{$b}->{rank}} keys %status_id_to_acc_list) {
	
	my $status_struct = $status_info{$status_id};
	my ($requires_update, $fails_incorporation, $status_descr) = ($status_struct->{requires_update},
								      $status_struct->{fails_incorporation},
								      $status_struct->{status_descr});
	my $accList_aref = $status_id_to_acc_list{$status_id};
	my $update_bgcolor = ($requires_update) ? "bgcolor=\"#00FF00\"" : "";
	my $incorp_bgcolor = ($fails_incorporation) ? "bgcolor=\"#FF0000\"" : "";
	my $accListing = "";
	foreach my $cdna_acc (@$accList_aref) {
	    $accListing .= "<a href=\"cdnaasmbl_report.cgi?cdna_acc=$cdna_acc&db=$db&compare_id=$compare_id_B\">$cdna_acc</a><br>\n";
	}
	
	$diff_text .= "<tr><td>$status_id</td><td>$status_descr</td><td $update_bgcolor>&nbsp;</td><td $incorp_bgcolor>&nbsp;</td><td>$accListing</td></tr>\n";
	
    }

    $diff_text .= "</table><p>";
}

