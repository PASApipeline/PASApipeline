#!/usr/bin/env perl

use Pasa_init;
use Pasa_conf;
use DBI;
use CGI;
use CGI::Pretty ":standard";
use CGI::Carp qw(fatalsToBrowser);
use Data::Dumper;
use Mysql_connect;
use Ath1_cdnas;
use strict;
use Gene_obj;
use Pasa_CGI;

my $cgi = new CGI;
print $cgi->header('text/html');

my $subcluster_id = $cgi->param('subcluster_id');
my $cdna_acc = $cgi->param('cdna_acc');
unless ($subcluster_id || $cdna_acc) {
    die "Need subcluster_id or cdna_acc\n";
}

my $db = $cgi->param('db') or die "Need db";
my $compare_id = $cgi->param('compare_id') or die "need compare_id";
our $SEE = ($cgi->param('SEE') || $cgi->param('DEBUG'));
my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $mysql_ro_user = &Pasa_conf::getParam("MYSQL_RO_USER");
my $mysql_ro_password = &Pasa_conf::getParam("MYSQL_RO_PASSWORD");
my ($dbproc) = &connect_to_db($mysql_server,$db,$mysql_ro_user,$mysql_ro_password);

my $css_common_text = &Pasa_CGI::get_common_CSS();

print $cgi->start_html(-title=>"cdna assembly report",
                       -head => style( { type => "text/css" }, $css_common_text ),
    );
my $annot_db = $cgi->param('annot_db');

## Custom URLs
my %link_templates;
&populate_link_templates();
my %global_link_substitutions = ( '__END5__' => '',   ## These are the supported attributes.
				  '__END3__' => '',   ##  Values are populated below
				  '__PASA_ACC__' => '',
				  '__CDNA_ACC__' => '',
				  '__GENE_ID__' => '',
				  '__MODEL_ID__' => '',
				  '__CONTIG_ID__'); 

my $form_counter = 0;


print <<_EOJS_;


<SCRIPT LANGUAGE=JAVASCRIPT 1.2>

function openWindow(url) {
    var props = "toolbar=yes,location=yes,directories=yes,status=yes,menubar=yes,width=" + window.innerWidth + ",height=" + window.innerHeight + ",scrollbars=yes,resizable=yes";
    win=window.open(url,"_blank",props);
}
</script>

_EOJS_

    ;


unless ($subcluster_id) {
    my $query = "select subcluster_id from subcluster_link sl where sl.cdna_acc = ? ";
    $subcluster_id = &very_first_result_sql($dbproc, $query, $cdna_acc);
    die "Couldn't get subcluster_id for cdna_acc" unless ($subcluster_id);
}

my $cluster_id = &Ath1_cdnas::get_cluster_id_via_subcluster_id($dbproc, $subcluster_id);
my $query = "select cluster_id, lend, rend, annotdb_asmbl_id from clusters where cluster_id = $cluster_id";
my $result = &first_result_sql($dbproc, $query);
my ($clusterid, $lend, $rend, $annotdb_asmbl_id) = @$result;

$global_link_substitutions{'__CONTIG_ID__'} = $annotdb_asmbl_id;
$global_link_substitutions{'__END5__'} = $lend;
$global_link_substitutions{'__END3__'} = $rend;


print "<h1>Report for cDNA subcluster: $subcluster_id</h1>";
print "<h3> of cluster: $cluster_id (annotdb_asmbl_id:$annotdb_asmbl_id coords:$lend-$rend)</h2>";

print &create_link_selector("contig_id", $annotdb_asmbl_id);

# Print image of the assembly.
print "<h2>Subcluster view.</h2>\n";
print "<img src =\"asmbl_image_generator.cgi?subcluster_id=$subcluster_id&db=$db&IMAGE_X_SIZE=750&DRAW_PANEL_SCALER=0.7&compare_id=$compare_id\" alt=asmbl_image >\n";


## describe the assemblies
print "<h2>Assembly description</h2>";

my %assemblies;
my $query = "select cdna_acc from subcluster_link where subcluster_id = $subcluster_id";
my @assemblies;
my @results = &do_sql_2D($dbproc, $query);
foreach my $result (@results) {
    my $cdna_acc = $result->[0];
    push (@assemblies, $cdna_acc);
}

## get the underlying cDNAs built into the assembly:
my %cdnas;
my %assembly_to_cdnas;
foreach my $assembly (@assemblies) {
    my $query = "select cdna_acc from asmbl_link where asmbl_acc = ?";
    my @results = &do_sql_2D ($dbproc, $query, $assembly);
    my @x;
    foreach my $result (@results) {
        my $cdna = $result->[0];
        $cdnas{$cdna}++;
        push (@x, $cdna);
    }
    $assembly_to_cdnas{$assembly} = [@x];
}

# print assembly listing
print "<table border = 3>\n";
print "<th>assembly</th><th>cdnas</th><th>annotations linked</th><th>status</th><th>comment</th></tr>\n";
foreach my $assembly (@assemblies) {
    
    $global_link_substitutions{'__PASA_ACC__'} = $assembly;
    
    print "<tr><td rowspan=2>$assembly</td>";

    print "<td>&nbsp;</td>"; # replace the long cdna list
    
    #print "<td><table>\n";
    #my @cdnas = @{$assembly_to_cdnas{$assembly}};
    #foreach my $cdna (@cdnas) {
    #    my $bg_color = ($cdnas{$cdna} > 1) ? "#FF0000" : "#00FF00";
    #    print "<tr><td bgcolor=\"$bg_color\">$cdna</td></tr>\n";
    #}
    #print "</table></td>\n";
    ## get the annotations
    my @annotations = Ath1_cdnas::get_annotations_linked_to_assembly($dbproc, $assembly, $compare_id);
    foreach my $annotation (@annotations) {
        my ($gene_id, $model_id) = @$annotation;
        $global_link_substitutions{__GENE_ID__} = $gene_id;
        $global_link_substitutions{__MODEL_ID__} = $model_id;
        
        $annotation = "<tr><td>gene: $gene_id " . &create_link_selector("gene_id", $gene_id);
		
        if ($model_id) {
            $annotation .=  "<br>&nbsp;model: $model_id" . &create_link_selector("model_id", $model_id);
            $annotation .= "<a href=\"gene_centric_image_generator.cgi?db=$db&compare_id=$compare_id&TU=$gene_id\" target=_new >AllGeneUpdates</a>";
        }
        $annotation .= "</td></tr>";
        
        
    }
    
    if (@annotations) {
        unshift (@annotations, "<table border=1>");
        push (@annotations, "</table>");
    }
    
    print "<td>@annotations</td>";
    my @status_info = &Ath1_cdnas::get_status_info($dbproc, $assembly, $compare_id);
    my $status = pop @status_info;
    my $update_id = $status->{update_id};
    my $protein = "";
    if ($update_id && (my $gene_obj = &Ath1_cdnas::get_after_gene_via_update_id($dbproc, $update_id))) {
        $protein = $gene_obj->{protein_seq};
    }
    
    my $update_info = Ath1_cdnas::get_annot_update_info($dbproc, $update_id);
    my $comment =  $status->{comment} . "<br>Applied to " . $update_info->{gene_id} . "," . $update_info->{model_id};
    if ($protein) {
        $protein =~ s/(\w{60})/$1\n/g;
        $comment .= "<br><pre>&gt;$assembly-based protein<br><code>$protein</code></pre>";
    }
    
    print "<td>" . $status->{status_id} . ". " . $status->{status_descr} . "</td><td>" . $comment . "</td></tr>\n";
    print "<tr><td colspan=4><img src=\"asmbl_image_generator.cgi?assembly_acc=$assembly&db=$db&compare_id=$compare_id&IMAGE_X_SIZE=700&DRAW_PANEL_SCALER=0.8\" alt=asmbl_image></td></tr>\n";
    print "<tr><td colspan=5 bgcolor=\"#000000\">&nbsp;</td></tr>\n";
}
print "</table>";


## Get the cDNA status info:
print "<h2>cDNA information</h2>\n";
print "<table border = 2>\n";
print "<tr><th>cDNA acc</th><th>alignment program</th><th>spliced orientation</th><th>Seq Length</th><th>Full-length status</th><th>Title</th></tr>\n";
foreach my $cDNA (keys %cdnas) {
    my $align_id = Ath1_cdnas::get_align_id_via_align_acc($dbproc, $cDNA);
    my $cdna_info = Ath1_cdnas::get_cdna_info_via_align_id($dbproc, $align_id);
    my $bgcolor = ($cdna_info->{is_fli}) ? "#0000FF" : "#FFFFFF";
    $global_link_substitutions{__CDNA_ACC__} = $cDNA;
    
    print "<tr><td rowspan=2>" . $cdna_info->{cdna_acc} . "</td><td>" . $cdna_info->{prog} . "</td><td>" . $cdna_info->{spliced_orient} . "</td><td>" . $cdna_info->{length} . "</td><td bgcolor = \"$bgcolor\">" . $cdna_info->{is_fli} . "</td><td>" . $cdna_info->{header} . "</td></tr>\n";
    print "<tr><th></th><th>avg_per_id: $cdna_info->{avg_per_id}</th><th>percent_aligned: $cdna_info->{percent_aligned}</th><td colspan=3>Alignment: $cdna_info->{alignment}</td></tr>\n";
}
print "</table>\n";

print $cgi->end_html();

$dbproc->disconnect;



####
sub create_link_selector {
    my ($url_var_name, $value) = @_;
    
    $form_counter++;
    
    my $form_name = "form$form_counter";
    my $select_name = "links$form_counter";
    
    my $links_aref = $link_templates{$url_var_name};
    
    my $selection_component = "";
    if (ref $links_aref) {
        $selection_component = "<form name=$form_name method=get " 
            . "onChange=\"if (document.$form_name.$select_name.selectedIndex != 0) { openWindow(document.$form_name.${select_name}\[document.$form_name.$select_name.selectedIndex].value); }\">\n"
            . "<select name=$select_name ><option>-$value-links-</option>\n";
        
        my @links = @$links_aref;
        foreach my $struct (@links) {
            my ($url_name, $url_template) = ($struct->{url_name}, $struct->{url_template});
            $url_template =~ s/__URLVAR__/$value/;
            
            foreach my $sub_field (keys %global_link_substitutions) {
                my $value = $global_link_substitutions{$sub_field};
                $url_template =~ s/$sub_field/$value/g;
            }
            
            $selection_component .= "<option value=\"$url_template\">$url_name</option>\n";
        }
        $selection_component .= "</select></form>\n";
        
    }
    return ($selection_component);
}


####
sub populate_link_templates {
    my $query = "select url_name, url_template, url_var_name from URL_templates";
    my @results = &do_sql_2D($dbproc, $query);
    
    foreach my $result (@results) {
        
        my ($url_name, $url_template, $url_var_name) = @$result;
        my $list_ref = $link_templates{$url_var_name};
        unless (ref $list_ref) {
            $list_ref = $link_templates{$url_var_name} = [];
        }
        my $url_struct = { url_name => $url_name,
                           url_template => $url_template };
        push (@$list_ref, $url_struct);
        
    }
    
}



