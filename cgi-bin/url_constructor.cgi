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
use CGI::Pretty ":standard";
use Pasa_CGI;



$|++;

my $cgi = new CGI;
print $cgi->header('text/html');

my $css_common_text = &Pasa_CGI::get_common_CSS();

print $cgi->start_html(-title=>"PASA URL Constructor",
                       -head => style( { type => "text/css" }, $css_common_text ),
                       );

my $db = $cgi->param('db') or die "Error, require parameter db\n";

my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $mysql_ro_user = &Pasa_conf::getParam("MYSQL_RO_USER");
my $mysql_ro_password = &Pasa_conf::getParam("MYSQL_RO_PASSWORD");
my ($dbproc) = &connect_to_db($mysql_server,$db,$mysql_ro_user,$mysql_ro_password);


my $default_url_template = "http://your.server.address/path/to/cgi_script?TOKEN=__URLVAR__";
my $new_url_template = $cgi->param("url_template");
my $new_url_varname = $cgi->param("url_var");
my $new_url_name = $cgi->param("url_name");

&process_removals();

if ($new_url_template && $new_url_template !~ /your\.server\.address/) {
    &add_new_url_template($new_url_varname, $new_url_name, $new_url_template);
}

print "<center>\n";
print "<h1>PASA Custom Hyperlinking Facility</h1>\n";
print "<form name=url_create_form action=url_constructor.cgi method=post >\n";
print "<input type=hidden name=db value=\"$db\" />\n";

print "<h2>Existing URL Templates</h2>\n";
&printCurrentURLs();

print "<h2>Add a new URL Template</h2>\n";
&printForm();
print "</form>\n";


exit(0);

sub printForm {

    ## provide the construction inputs:
    my $query = "select url_var_name from URL_var_names";
    my @results = &do_sql_2D($dbproc, $query);
    my @var_names;
    foreach my $result (@results) {
	my ($var_name) = @$result;
	push (@var_names, $var_name);
    }
    
    print "<table border=1 >  " 
	. "<tr><th>url&nbsp;name</th><th>url&nbsp;var&nbsp;name</th><th>url&nbsp;template</th></tr>\n"
	
	. "<tr><td><input type=text name=url_name size=20 />"
	
	. "<td><select name=url_var >";
    
    foreach my $var_name (@var_names) {
	print "<option value=$var_name>$var_name</>";
    }
    
    print "</select></td>"
	
	
	
	. "<td><input type=text name=url_template size=75 value=\"$default_url_template\"></td></tr></table>\n";
    
    print "<br><br><input type=submit value=\"Process Request\">\n";
    
    print "</CENTER>";
    print "<br><br><p>Custom URLs are constructed to provide linking capabilities from PASA components to external resources.  For example, given a gene_id that is linked to a PASA alignment assembly, a link can be constructed to an external gene report page (ie. Manatee).  The <b>url var name</b> identifies the PASA feature for which customized links can be provided.  The templated URL parameter &quot;__URLVAR__&quot; is replaced with the value of the <b>url var name</b> on a corresponding PASA report page.</p>\n";
    print "<p>Additional attributes that will be substituted in all URLs, when values are known, include:</p>"
	. "<ul>"
	. "<li>__END5__</li>"
	. "<li>__END3__</li>"
	. "<li>__PASA_ACC__</li>"
	. "<li>__CDNA_ACC__</li>"
	. "<li>__GENE_ID__</li>"
	. "<li>__MODEL_ID__</li>"
	. "<li>__CONTIG_ID__</li>"
	. "</ul>";
    
    
    print $cgi->end_html();
    
}


####
sub add_new_url_template {
    my ($new_url_varname, $new_url_name, $new_url_template) = @_;

    unless ($new_url_name) {
	&error("Sorry, no url name was specified, couldn't create new url template for:<br>$new_url_template");
	return;
    }
    unless ($new_url_template =~ /__URLVAR__/) {
	&error("Sorry, your url template:<br>$new_url_template<br>is missing the __URLVAR__ token.");
	return;
    }

    my $query = "insert URL_templates (url_name, url_template, url_var_name) values (?,?,?)";
    &RunMod($dbproc, $query, $new_url_name, $new_url_template, $new_url_varname);
    
}

####
sub error {
    my $error_text = shift;
    print "<font color=\"#ff0000\">$error_text</font>";
}

####
sub printCurrentURLs {
    my $query = "select url_name, url_template, url_var_name from URL_templates";
    my @results = &do_sql_2D($dbproc, $query);
    
    print "<table border=1 width=600>";
    print "<tr><th>delete</th><th>url&nbsp;name</th><th>url&nbsp;var&nbsp;name</th><th>url&nbsp;template</th></tr>\n";
    if (@results) {
	
	foreach my $result (@results) {
	    my ($url_name, $url_template, $url_var_name) = @$result;
	    print "<tr><td><input type=checkbox name=delete_$url_name /></td><td>$url_name</td><td>$url_var_name</td><td>" . CGI::escapeHTML($url_template) . "</td></tr>\n";
	}
    } else {
	print "<tr><td colspan=4 align=center><i> no existing url templates available</i></td></tr>\n";
    }
    print "</table><br><br>\n";
}


####
sub process_removals() {
    my %vars = $cgi->Vars();
    foreach my $key (%vars) {
	if ($key =~ /^delete_(\S+)/) {
	    my $remove_url_name = $1;
	    my $query = "delete from URL_templates where url_name = ?";
	    &RunMod($dbproc, $query, $remove_url_name);
	    &add_info("-removed url ($remove_url_name)<br>");
	}
    }
}


####
sub add_info {
    my $info = shift;
    print "<font color=\"#0000ff\">$info</font>\n";
}


    
