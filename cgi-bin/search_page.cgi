#!/usr/bin/env perl

use Pasa_init;
use Pasa_conf;
use strict;
use DBI;
use CGI;
use CGI::Pretty ":standard"; 
use CGI::Carp qw(fatalsToBrowser);
use Data::Dumper;
use Mysql_connect;
use Ath1_cdnas;
use Pasa_CGI; 

$|++;

our $DEBUG = 0;
my $cgi = new CGI;
print $cgi->header('text/html');
my $db = $cgi->param('db');
our $SEE = $cgi->param('SEE') || $cgi->param('DEBUG');
my $cdna_acc = $cgi->param('cdna_acc');
my $assembly_acc = $cgi->param('assembly_acc');
my $annot_id = $cgi->param('annot_id');
my $compare_id = $cgi->param('compare_id') || 1;
my $sub_match = $cgi->param('sub_match') || 0;
unless ($db) {
    die "Must set the db parameter.\n";
}

my $css_text = &Pasa_CGI::get_CSS($0);
print $cgi->start_html(-title=>"search page", 
                       -head => style( { type => "text/css" }, $css_text ),      
    );


my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $mysql_ro_user = &Pasa_conf::getParam("MYSQL_RO_USER");
my $mysql_ro_password = &Pasa_conf::getParam("MYSQL_RO_PASSWORD");
my ($dbproc) = &connect_to_db($mysql_server,$db,$mysql_ro_user,$mysql_ro_password);


if ($cdna_acc || $assembly_acc || $annot_id) {

    print "<h1>Search Results" . ($cdna_acc, $assembly_acc, $annot_id) . "</h1>\n";
    
    if ($cdna_acc) {
        &print_cdna_search_results($cdna_acc);
    } elsif ($assembly_acc) {
        &print_assembly_search_results($assembly_acc);
    } elsif ($annot_id) {
        &print_annot_id_search_results($annot_id);
    }
}

&print_input_form();

print $cgi->end_html();

$dbproc->disconnect;

exit(0);


####
sub print_input_form() {

    my $query = "select compare_id from annotation_compare order by compare_id asc";
    my @compare_ids;
    foreach my $result (&do_sql_2D($dbproc,$query)) {
        my ($id) = @$result;
        push (@compare_ids, $id);
    }
    
    
    my $compare_id_selector = "<select name='compare_id'>";
    {  # build the selector:
        if (@compare_ids) {
            foreach my $id (@compare_ids) {
                $compare_id_selector .= "<option value=\"$id\">$id</option>\n";
            }
        }
        else {
            ## nothing to choose.  Just set it to 1 for now.
            $compare_id_selector .= "<option value=\"1\">1</option>\n";
        }
        $compare_id_selector .= "</select>\n";
    }
    
    print "<h1>Search $db for PASA Alignment Assemblies</h1>\n";
    print <<__EOform;
    
    <FORM id="alignmentAssemblySearchForm" ACTION="search_page.cgi" METHOD="get">
        <input type=hidden name=db value="$db"/>
        <p>Alignment assembly accession:
        <input type=text name="assembly_acc"></p>    

    <p>cDNA accession:
    <input type=text name="cdna_acc"/>
<input type=checkbox name=sub_match value=1> Match accession substring<br><i> otherwise requires full-match (ie. 'gi|19871417|gb|AV829357.1|AV829357')</i></p>    
    
    <p>Annotation ID:
    <input type=text name="annot_id"/>
    Compare_id: $compare_id_selector</p>

    <input type=submit value="Search"/>
	</FORM>
	
__EOform

    ;
    
};


## Search cDNA accession:
sub print_cdna_search_results {
    my $cdna_acc = shift;
    print "<h2>search results for cDNA accesssion: $cdna_acc</h2>\n";
    my $query = "select asmbl_acc from asmbl_link where cdna_acc = \"$cdna_acc\"";
    if ($sub_match) {
	$query = "select asmbl_acc, cdna_acc from asmbl_link where cdna_acc like \"\%$cdna_acc\%\"";
    }
    my @results = &do_sql_2D($dbproc, $query);
    my @asmbl_accs;
    if (@results) {
        foreach my $result (@results) {
            my ($asmbl_acc, $cdna_acc) = @$result;
            print "<p class='searchResults'>Found cDNA accession <span class='searchToken'>$cdna_acc</span> in assembly $asmbl_acc</p>\n";
            &print_assembly_search_results($asmbl_acc);
        }
    } else {
        print "<p class='noSearchResults'>Sorry, no matches for cDNA accession: <span class='searchToken'>$cdna_acc</span></p>\n";
    }
    
}


####
sub print_assembly_search_results {
    my $assembly_acc = shift;
    my $query = "select subcluster_id from subcluster_link where cdna_acc = ?";
    my $result = &very_first_result_sql($dbproc, $query, $assembly_acc);
    if (my $subcluster_id = $result) {
        print "<p class='searchResults'>Found assembly <span class='searchToken'>$assembly_acc </span>in subcluster: <a href=\"cdnaasmbl_report.cgi?subcluster_id=${subcluster_id}&db=${db}&compare_id=${compare_id}\">$subcluster_id</a></p>\n";
    }
}




## Search for annotated genes
sub print_annot_id_search_results {
    my $annot_id = shift;
    print "<h2>Search results for annotation ID: <span class='searchToken'>$annot_id</span></h2>\n";
    my $query = "select distinct cdna_acc from annotation_link where compare_id = $compare_id and (model_id = \"$annot_id\" or gene_id = \"$annot_id\")";
    my @cdna_accs;
    my @results = &do_sql_2D($dbproc, $query);
    foreach my $result (@results) {
        my ($cdna_acc) = @$result;
        push (@cdna_accs, $cdna_acc);
    }
    if (@cdna_accs) {
        foreach my $cdna_acc (@cdna_accs) {
            &print_assembly_search_results($cdna_acc);
        }
    } else {
        print "<p class='noSearchResults'>Sorry, no search results for annotation_id: <span class='searchToken'>$annot_id</span></p>\n";
    }
}

