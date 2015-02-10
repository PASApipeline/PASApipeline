#!/usr/bin/env perl

use Pasa_init;
use Pasa_conf;
use strict;
use DBI;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Data::Dumper;
use Mysql_connect;
use Ath1_cdnas;

$|++;

our $DEBUG = 0;
my $cgi = new CGI;
print $cgi->header('text/html');
my $db = $cgi->param('db');
our $SEE = $cgi->param('SEE') || $cgi->param('DEBUG');

unless ($db) {
    die "Must set the db parameter.\n";
}

print $cgi->start_html(-title=>"search for alternative splicing variants");
my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $mysql_ro_user = &Pasa_conf::getParam("MYSQL_RO_USER");
my $mysql_ro_password = &Pasa_conf::getParam("MYSQL_RO_PASSWORD");
my ($dbproc) = &connect_to_db($mysql_server,$db,$mysql_ro_user,$mysql_ro_password);

if ($cgi->param("SEARCH")) {
    &get_search_results();
}
else {
    &print_search_form();
}

print $cgi->end_html();

$dbproc->disconnect;

exit(0);



####
sub get_search_results {
    my %params = $cgi->Vars();
    
    my $num_subfeatures = $params{numSubVariations};
    
    ## unpack the variation types
    my @variation_types = split (/\0/, $params{spliceVariationTypes});
    
    my @cdna_accs;


    if (@variation_types) {
        print "<h1>Search results for: " . join (", ", @variation_types) . "</h1>\n";
        ## get all cdna accessions that have all types
        my %cdna_to_types;
        my $type_list = join ('","', @variation_types);
        my $query = "select cdna_acc, type from splice_variation where type in (\"$type_list\")";
        my @results = &do_sql_2D($dbproc, $query);
        foreach my $result (@results) {
            my ($cdna_acc, $type) = @$result;
            $cdna_to_types{$cdna_acc} .= " $type ";
        }
        
        ## check to see that each has all the requested types
        foreach my $cdna (keys %cdna_to_types) {
            my $type_listing = $cdna_to_types{$cdna};
            my $got_all = 1;
            foreach my $type (@variation_types) {
                unless ($type_listing =~ /\b$type\b/) {
                    $got_all = 0;
                    last;
                }
            }
            if ($got_all) {
                push (@cdna_accs, $cdna);
            }
        }
        
        if (scalar (@variation_types) == 1 && $num_subfeatures ne "any") {
            ## check for number of splice variations:
            print "<h2>Restricted to those with number of event features = $num_subfeatures</h2>\n";
            my @cdna_accs_copy = @cdna_accs;
            @cdna_accs = (); #re init
            
            foreach my $cdna (@cdna_accs_copy) {
                my $query = "select num_subfeatures_included from splice_variation where type = \"$variation_types[0]\" and cdna_acc = \"$cdna\" order by num_subfeatures_included\n";
                my @results = &do_sql_2D($dbproc, $query);
                foreach my $result (@results) {
                    my ($num_sub) = @$result;
                    if ($num_sub eq $num_subfeatures) {
                        push (@cdna_accs, $cdna);
                        last;
                    }
                }
            }
        }
    }
    elsif (my $cdna = $params{cdna_acc}) {
        print "<h1>Search results for accession $cdna</h1>\n";
        my $query = "select count(*) from splice_variation where cdna_acc = \"$cdna\"";
        my $count = &very_first_result_sql($dbproc, $query);
        if ($count) {
            push (@cdna_accs, $cdna);
        }
    }
    
    if (@cdna_accs) {
        &print_cdna_listing(@cdna_accs);
    }
    else {
        print "<font color=\"FF0000\">Sorry, no entries found.</font>\n";
    }
}

####
sub print_search_form {
    
    my $query = "select distinct type from splice_variation";
    my @results = &do_sql_2D($dbproc, $query);
    my @types;
    foreach my $result (@results) {
        my $type = $result->[0];
        push (@types, $type);
    }

    print "<h1>Search for Alternative Splicing Isoforms</h2>\n";
    print "<form>\n";
    print "<ul>\n";
    
    ## build type selection list
    print "<li>Splicing Variation: ";
    print "<select name=spliceVariationTypes size=9 multiple>\n";
    foreach my $type (sort @types) {
        print "<option value=$type>$type</option>\n";
    }
    print "</select>\n";

    ## number of subfeatures
    my @distinct_subfeature_counts;
    my $query = "select distinct num_subfeatures_included from splice_variation order by num_subfeatures_included";
    @results = &do_sql_2D($dbproc, $query);
    foreach my $result (@results) {
        my $num_subfeats = $result->[0];
        push (@distinct_subfeature_counts, $num_subfeats);
    }
    print "<li>Count of subfeatures: ";
    print "<select name=numSubVariations ><option value=any selected>any</option>\n";
    foreach my $subfeat_count (@distinct_subfeature_counts) {
        print "<option value=$subfeat_count>$subfeat_count</option>\n";
    }
    print "</select> ";
    print "<i>if multiple splicing variation types are chosen, this setting is ignored.</i>\n";
    
    print "<li>PASA assembly accession: <input type=text name=cdna_acc /> ";
        

    print "</ul>\n";
    
    print "<input type=submit value=\"Find Entries\">\n";
    print "<input type=reset value=\"Reset\">\n";
    print "<input type=hidden name=db value=$db >\n";
    print "<input type=hidden name=SEARCH value=1 >";
    
    print "</form>\n";

}

####
sub print_cdna_listing {
    my @cdnas = @_;
    print "<table border=1>\n";
    foreach my $cdna (@cdnas) {
        print "<tr><td><a href=\"assembly_alt_splice_info.cgi?db=$db&cdna_acc=$cdna\">$cdna</a></td></tr>\n";
    }
    print "</table>\n";
}
