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
use Pasa_CGI;

$|++;

my $cgi = new CGI;
print $cgi->header('text/html');
my $db = $cgi->param('db');
our $SEE = $cgi->param('SEE') || $cgi->param('DEBUG');
our $DB_SEE = $SEE;


my $no_cache = $cgi->param('NOCACHE');

unless ($db) {
    die "Must set the db parameter.\n";
}

my $css_common_text = &Pasa_CGI::get_common_CSS();
$css_common_text .= &Pasa_CGI::get_page_specific_CSS($0);


unless ($ENV{WEBSERVER_TMP}) { $ENV{WEBSERVER_TMP} = "/tmp";}

if ($SEE) {
    print "<PRE>\n";
}


my $cache_file = "$ENV{WEBSERVER_TMP}/$db.status.html";
my $time1 = time();
print $cgi->start_html(-title=>"PASA-2 status report for $db",
                       -head => style( { type => "text/css" }, $css_common_text ),
    );

my $javascript_code = <<_EOJS_;

<script type="text/javascript">
<!--
    function toggleVisibility(divID, textID) {
        

      
       var e = document.getElementById(divID);
       var t = document.getElementById(textID); 
       
       if(e.style.display == 'block') {
          e.style.display = 'none';
          t.innerHTML = "+";
       }
       else {
          e.style.display = 'block';
          t.innerHTML = "-";
       }
     
    }
//-->
</script>

_EOJS_

        ;


print $javascript_code;



## use cached data whenever possible.  It takes some serious time to generate this file for large data sets.
if (-s $cache_file && ! $no_cache) {
    print "<p id='refreshLink'>Click <a href=\"status_report.cgi?db=$db&NOCACHE=1\">here</a> to freshly reload the data only as needed. This page was cached.</p>\n";
    system "cat $cache_file";
    print $cgi->end_html();
    exit(0);
}

## write the file

open (CACHE, ">$cache_file") or die "Error, cannot write to $cache_file";

my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $mysql_ro_user = &Pasa_conf::getParam("MYSQL_RO_USER");
my $mysql_ro_password = &Pasa_conf::getParam("MYSQL_RO_PASSWORD");
    
my ($dbproc) = &connect_to_db($mysql_server,$db,$mysql_ro_user,$mysql_ro_password);


## Simple way to keep users from using this script 
my $query = "select count(*) from status";
my $status_count = &very_first_result_sql($dbproc, $query);
if ($status_count < 44) {
    die "Sorry, this PASA database is incompatible with the latest version of PASA. ";
}


my %requires_update;
my $query = "select status_id, requires_update from status";
my @results = &do_sql_2D($dbproc, $query);
foreach my $result (@results) {
    my ($status_id, $requires_update) = @$result;
    $requires_update{$status_id} = $requires_update;
}





my $text = "<h1>PASA Database $db Report</h1>\n";
# Get the number of cdnas:
my $query = "select count(*) from cdna_info where is_assembly = 0";
my $total_cdnas = &very_first_result_sql($dbproc, $query);

# Get the number of fli-cdnas
my $query = "select count(*) from cdna_info where is_assembly = 0 and is_fli = 1";
my $total_fli = &very_first_result_sql($dbproc, $query);

# Get the number of non-fli cdnas:
my $query = "select count(*) from cdna_info where is_assembly = 0 and is_fli = 0";
my $total_non_fli = &very_first_result_sql($dbproc, $query);

my @prog_and_valid_aligns;

## number cdnas with any alignment (number mapped)
my $query = "select count(distinct ci.cdna_acc) from cdna_info ci, align_link al where ci.id = al.cdna_info_id and ci.is_assembly = 0";
my $number_mapped_transcripts = &very_first_result_sql($dbproc, $query);


## valid alignments:
my $query = "select al.prog, count(distinct ci.cdna_acc) from cdna_info ci, align_link al where al.validate = 1 and al.cdna_info_id = ci.id and ci.is_assembly = 0 group by prog order by 2 desc";
my @results = &do_sql_2D($dbproc, $query);
foreach my $result (@results) {
    my ($prog, $count) = @$result;
    push (@prog_and_valid_aligns, [$prog, $count]);
}


# Get the total number of validating cdnas:
my $query = "select count(distinct ci.cdna_acc) from cdna_info ci, align_link al where al.validate = 1 and al.cdna_info_id = ci.id and ci.is_assembly = 0";
my $total_valid = &very_first_result_sql($dbproc, $query);

## get total number of validating FL-cDNA alignments
my $query = "select count(distinct ci.cdna_acc) from cdna_info ci, align_link al where al.validate = 1 and al.cdna_info_id = ci.id and ci.is_assembly = 0 and ci.is_fli = 1";
my $total_valid_FL_cdnas = &very_first_result_sql($dbproc, $query);

## get total number of validating non-FL alignments
my $query = "select count(distinct ci.cdna_acc) from cdna_info ci, align_link al where al.validate = 1 and al.cdna_info_id = ci.id and ci.is_assembly = 0 and ci.is_fli = 0";
my $total_valid_ESTs = &very_first_result_sql($dbproc, $query);


# Get the total number of clusters:
my $query = "select count(distinct cluster_id) from align_link where validate = 1";
my $num_clusters = &very_first_result_sql($dbproc, $query);

# Get the total number of subclusters:
my $query = "select max(subcluster_id) from subclusters";
my $num_subclusters = &very_first_result_sql($dbproc, $query);

# Get the total number of assemblies:
my $query = "select count(*) from cdna_info where is_assembly = 1";
my $num_assemblies = &very_first_result_sql($dbproc, $query);

# Get the total number of assemblies containing fli-cdnas:
my $query = "select count(*) from cdna_info where is_assembly = 1 and is_fli = 1";
my $num_fli_assemblies = &very_first_result_sql($dbproc, $query);

# Get the total number of assemblies containing non-fli cdnas:
my $query = "select count(*) from cdna_info where is_assembly = 1 and is_fli = 0";
my $num_non_fli_assemblies = &very_first_result_sql($dbproc, $query);

$text .= "<table id=\"pasaAlignmentAssemblyStatsTable\">"
    . "<tr><th>Transcripts or Assemblies</th><th>Count</th></tr>\n"
    . "<tr><td>Total transcript seqs</td><td>$total_cdnas</td></tr>\n"
    . "<tr><td>Fli cDNAs</td><td>$total_fli</td></tr>"
    . "<tr><td>partial cDNAs (ESTs)</td><td>$total_non_fli</td></tr>"
    . "<tr><td>Number transcripts with any alignment</td><td>$number_mapped_transcripts</td></tr>\n";

foreach my $prog_and_count (@prog_and_valid_aligns) {
    my ($prog, $count) = @$prog_and_count;
    $text .= "<tr><td>Valid <b>$prog alignments</b></td><td>$count</td></tr>\n";
}

$text .=   "<tr><td>Total Valid alignments</td><td>$total_valid</td></tr>" 
    . "<tr><td>Valid FL-cDNA alignments</td><td>$total_valid_FL_cdnas</td></tr>\n"
    . "<tr><td>Valid EST alignments</td><td>$total_valid_ESTs</td></tr>\n"
    
    . "<tr><td>Number of assemblies</td><td>$num_assemblies</td></tr>"
    . "<tr><td>Number of subclusters (genes)</td><td>$num_subclusters</td></tr>"
    . "<tr><td>Number of fli-containing assemblies</td><td>$num_fli_assemblies</td></tr>"
    . "<tr><td>Number of non-fli-containing assemblies</td><td>$num_non_fli_assemblies</td></tr></table>\n";


##################################
## Links to other pasa utilities

$text .= "<h2>PASA resources:</h2>\n";
$text .= "<ul id=\"pasaLinksList\">\n";

$text .= "<li><a href=\"describe_alignment_assemblies.cgi?db=$db\" target=_blank >Describe alignment assemblies</a>\n";

$text .= "<li><a href=\"describe_subclusters.cgi?db=$db\" target=_blank >Describe subclusters of assemblies</a>\n";

$text .= "<li><a href=\"retrieve_assembly_sequences.cgi?db=$db\" target=_blank >Retrieve alignment assembly tentative cDNA sequences</a>\n";


$text .= "<li><a href=\"search_page.cgi?db=$db\" target=\"search_$db\">Search the database.</a>\n";

$text .= "<li><a href=\"url_constructor.cgi?db=$db\" target=\"url_constructor_$db\">Construct customized URLs linked from PASA assembly report pages.</a>";

$text .= "<li><a href=\"alt_splice_report.cgi?db=$db\" target= \"alternative_splicing\">Alternative Splicing Report</a>";
$text .= "</ul>\n";


######################################################################################################################
########################## Results for Each Comparison ###############################################################
######################################################################################################################


$text .= "<h1>Status Report for cDNA Assembly Annotation Comparison</h1>\n";

my $query = "select compare_id, annotation_version, date from annotation_compare order by compare_id desc";
my @results = &do_sql_2D($dbproc, $query);
my $table_counter = 0;
foreach my $result (@results) {
    my ($compare_id, $annotation_version, $date) = @$result;
    
	$table_counter++;
	
	my $view_token = "+"; #($table_counter == 1) ? '-' : '+';
	$text .= "<div id=\"compare_$compare_id\" class=annotComparison>\n"; 
    $text .= "<h2><a href='#' onclick=\"javascript:toggleVisibility(\'compare_${compare_id}_contents\', \'TextToggle_${compare_id}\');\">[<span id=\'TextToggle_${compare_id}\'>$view_token</span>]</a>";
	$text .= " Comparison # $compare_id, Annotation Version: $annotation_version, date: $date</h2>\n";
    
	my $display_type = "none"; #($table_counter == 1) ? "block" : "none";
    $text .= "<div id=\"compare_${compare_id}_contents\" style=\"display: $display_type;\" >\n";
	
    $text .= "<i>some FL-assemblies may have been deprecated to EST status</i>\n";
    
    my %status_summary;
    
    my $query = "select s.status_id, s.requires_update, s.fails_incorporation, s.status_descr from status s order by s.rank";
    
    my @results = &do_sql_2D($dbproc, $query);
    
    foreach my $result_aref (@results) {
        my ($status_id, $requires_update, $fails_incorporation, $status_descr) = @$result_aref;
        
        
        my $query = "select count(*) from status_link sl where sl.status_id = $status_id and sl.compare_id = $compare_id";
        my $count = &very_first_result_sql($dbproc, $query);
        
        
        
        $status_summary{$status_id} = { count => $count,
                                        requires_update => $requires_update,
                                        fails_incorporation => $fails_incorporation,
                                        descr => $status_descr };
        
    }
    
    
    $text .= &assembly_status_table($compare_id, \%status_summary);
    
    $text .= "<br>";
    
    $text .= &assembly_incorporation_summary($dbproc, $compare_id);
    
    
    $text .= "<p><a href=\"passAndFail_GeneLists.cgi?db=$db&compare_id=$compare_id\" target=_new >Retrieve Lists of Genes Corresponding to Pass and Fail Categories.</a></p>";
    
    
    ############## Section Annotation Updates Summary  #########################################################
    
    ## get annotation updates.
    $text .= "<h2>Proposed Annotation Updates <i><font size=-1>[counts as unique models, not genes]</font></i></h2>";
    
	$text .= "<table><th>status_id</th><th>Description</th><th>Num Gene Model Updates</th><th>Num Alt Splice isoforms to Add</th><th>Num Novel Genes to Add</th></tr>\n";
    
    foreach my  $result_aref (@results) {
        my ($status_id, $requires_update, $fails_incorporation, $status_descr, $count) = @$result_aref;
        unless ($requires_update{$status_id}) { next;}
        my $query = "select count(distinct a.update_id) from status_link sl, annotation_updates a where sl.status_id = $status_id and sl.compare_id = $compare_id and sl.annot_update_id = a.update_id and a.compare_id = $compare_id and a.is_valid = 1 and a.alt_splice_flag = 0 and a.is_novel_flag = 0 and a.model_id is not NULL";
        my $total_gene_updates = &very_first_result_sql($dbproc, $query);
        
        my $query = "select count(distinct a.update_id) from status_link sl, annotation_updates a where sl.status_id = $status_id and sl.compare_id = $compare_id and sl.annot_update_id = a.update_id and a.compare_id = $compare_id and a.is_valid = 1 and a.alt_splice_flag = 1";
        my $total_alt_splice = &very_first_result_sql($dbproc, $query);
        
        my $query = "select count(distinct a.update_id) from status_link sl, annotation_updates a where sl.status_id = $status_id and sl.compare_id = $compare_id and sl.annot_update_id = a.update_id and a.compare_id = $compare_id and a.is_valid = 1 and a.alt_splice_flag = 0 and a.is_novel_flag = 1";
        my $total_new_genes = &very_first_result_sql($dbproc, $query);
        
        $text .= "<tr><td>$status_id</td><td><a href=\"gene_update_report.cgi?db=$db&status_id=$status_id&compare_id=$compare_id\">$status_descr</a></td><td>$total_gene_updates</td><td>$total_alt_splice</td><td>$total_new_genes</td></tr>\n";
        
    }
    
    my $query = "select count(distinct a.update_id) from status_link sl, annotation_updates a, status s where sl.compare_id = $compare_id and sl.status_id = s.status_id and s.requires_update = 1 and sl.annot_update_id = a.update_id and a.compare_id = $compare_id and a.is_valid = 1 and a.alt_splice_flag = 0 and a.is_novel_flag = 0 and a.model_id is not NULL";
    my $update_counter = &very_first_result_sql($dbproc, $query);
    
    my $query = "select count(distinct a.update_id) from status_link sl, annotation_updates a, status s where sl.compare_id = $compare_id and sl.status_id = s.status_id and s.requires_update = 1 and sl.annot_update_id = a.update_id and a.compare_id = $compare_id and a.is_valid = 1 and a.alt_splice_flag = 1";
    my $alt_splice_counter = &very_first_result_sql($dbproc, $query);
    
    my $query = "select count(distinct a.update_id) from status_link sl, annotation_updates a, status s where sl.compare_id = $compare_id and sl.status_id = s.status_id and s.requires_update = 1 and sl.annot_update_id = a.update_id and a.compare_id = $compare_id and a.is_valid = 1 and a.alt_splice_flag = 0 and a.is_novel_flag = 1";
    my $new_gene_counter = &very_first_result_sql($dbproc, $query);
    
    $text .= "<tr><th>&nbsp;</th><th>Totals <i><font=-2>(some models in multiple classes)</font></i></th><th>$update_counter</th><th>$alt_splice_counter</th><th>$new_gene_counter</th></tr></table>\n";


	# determine the number of protein sequences that have changed:
	my $query = "select count(distinct a.update_id) "
		. "from status_link sl, annotation_updates a, status s "
		. "where sl.compare_id = $compare_id "
		. "and sl.status_id = s.status_id "
		. "and s.requires_update = 1 "
		. "and sl.annot_update_id = a.update_id "
		. "and a.compare_id = $compare_id "
		. "and a.is_valid = 1 "
		. "and a.alt_splice_flag = 0 "
		. "and s.status_id not in (4,9, 10, 24, 33, 44, 13, 17, 25, 27, 32) "; # ignore alt splicing, novel, and UTR updates.

	my $num_prots_changed = &very_first_result_sql($dbproc, $query);
	
	$text .= "<p><b>Number of annotated proteins changed: $num_prots_changed</b>\n";
	

    $text .= "<p><a href=\"summarize_updates.cgi?db=$db&compare_id=$compare_id\">Summarize Gene Structure Updates to be Performed.</a></p>\n";
    
    
    #last;  ## For Debugging
    
    $text .= "</div>\n</div>\n<!-- end of annotation comparison $compare_id -->\n<p>";
    
}	


$dbproc->disconnect;
print CACHE $text;
close CACHE;

if ($SEE) {
    print "</PRE>\n";
}
 

print $text;

my $time2 = time();
my $time_to_load_page = $time2-$time1;
print "<p>Page originally generated in $time_to_load_page seconds.</p>\n";

print $cgi->end_html();

exit(0);


####
sub assembly_status_table {
    my ($compare_id, $status_summary_href) = @_;
    my %status_summary = %$status_summary_href;
    my $table_html = <<_EOTABLE;
    
    <table class="alignmentAssemblyClassificationsTable">
	<tr><th colspan=5>Annotation Classification for Alignment Assemblies</th></tr> 
    <tr><th></th><th colspan=2>FL-assemblies</th><th colspan=2>EST-assemblies</th></tr>
	<tr><th></th><th>PASS</th><th>fail</th><th>PASS</th><th>fail</th></tr>
	<tr><td>Incorporated</td><td>__STATUS_3__</td><td></td><td>__STATUS_12__</td><td></td></tr>
	<tr><td>UTR addition</td><td>__STATUS_4__</td><td></td><td>__STATUS_13__</td><td></td></tr>
	<tr><td>Gene extension</td><td>__STATUS_6__</td><td>__STATUS_7__</td><td>__STATUS_14__</td><td>__STATUS_15__</td></tr>
	<tr><td>Internal gene structure rearrangement</td><td></td><td>__STATUS_23__</td><td></td><td>__STATUS_18__</td></tr>
	<tr><td>__TAB__-passes homology tests</td><td>__STATUS_8__</td><td></td><td>__STATUS_16__</td><td></td></tr>
	<tr><td>__TAB__-fails homology, passes ORF span</td><td>__STATUS_26__</td><td></td><td>__STATUS_31__</td><td></td></tr>
	<tr><td>Gene Merging</td><td>__STATUS_29__</td><td>__STATUS_30__</td><td>__STATUS_36__</td><td>__STATUS_35__</td></tr>
	<tr><td>Gene Splitting</td><td>__STATUS_40__</td><td>__STATUS_41__</td><td></td><td></td></tr>
	<tr><td>Alt Splicing Isoform</td><td></td><td>__STATUS_21__</td><td></td><td></td></tr>
	<tr><td>__TAB__-passes homology test</td><td>__STATUS_9__</td><td></td>
    <td><table width=100% ><tr><td>__STATUS_25__</td></tr><tr><td>__STATUS_24__</td></tr></table>
    </td><td></td></tr>
	<tr><td>__TAB__-fails homology, passes ORF span</td><td>__STATUS_33__</td><td></td><td>__STATUS_32__</td><td></td></tr>
	<tr><td>New Gene</td><td>__STATUS_10__</td><td>__STATUS_22__</td><td></td><td>__STATUS_19__</td></tr>
	<tr><td>__TAB__Alt splice of new gene</td><td>__STATUS_44__</td><td>__STATUS_45__</td><td>__STATUS_27__</td><td>__STATUS_28__</td></tr>
	<tr><td colspan=5>&nbsp;</td></tr>
	<tr><td>FL-assembly fails gene requirements</td><td></td><td>__STATUS_1__</td><td></td><td></td></tr>
	<tr><td>Antisense</td><td></td><td>__STATUS_42__</td><td></td><td>__STATUS_43__</td></tr>
	<tr><td>Single-exon EST-assembly incompatible</td><td></td><td></td><td></td><td>__STATUS_34__</td></tr>
	<tr><td colspan=5>&nbsp;</td></tr>
	<tr><td>delayed incorporation due to gene merging</td><td></td><td>__STATUS_37__</td><td></td><td>__STATUS_38__</td></tr>
	<tr><td>delayed incorporation due to gene splitting</td><td></td><td>__STATUS_39__</td><td></td><td></td></tr>
	<tr><td colspan=5>&nbsp;</td></tr>
	<tr><td>Total</td><td colspan=4>__TOTAL__</td></tr>
    </table>
	

_EOTABLE
    
    
    ;

my $tab = "&nbsp;" x 5;
$table_html =~ s/__TAB__/$tab/g;

my $total = 0;
my $errors;
foreach my $status_id (sort keys %status_summary) {
	
	my $href = $status_summary{$status_id};
	
	my ($count, $requires_update, $fails_incorporation, $status_descr) = ($href->{count},
                                                                          $href->{requires_update},
                                                                          $href->{fails_incorporation},
                                                                          $href->{descr});
	$total += $count;
	
	my $bgcolor = "";
	if ($requires_update) {
	    $bgcolor = "bgcolor=\"#00FF00\"";
	} elsif ($fails_incorporation) {
	    $bgcolor = "bgcolor=\"#ff9f9f\"";
	}
	
	my $key = "<td>__STATUS_${status_id}__</td>";
    
    my $html_replacement = "<td $bgcolor><a href=\"status_to_asmbl_listing.cgi?status_id=$status_id&db=$db&compare_id=$compare_id\" target=_new  title=\"$status_descr\" ><b>$count</b></a></td>";
    
	my $substitution  = $table_html =~ s/$key/$html_replacement/;
	if ($count && !$substitution) {
	    $errors .= "<pre>There was no sub for $key, $count\n</pre>\n";
	}
}


$table_html =~ s/__TOTAL__/<b>$total<\/b>/;

$table_html .= $errors;

return ($table_html);

}


####
sub assembly_incorporation_summary {
    my ($dbproc, $compare_id) = @_;
    ## build model_to_TU_link;
    my %model_to_TU;
    my $query = "select gene_id, model_id from annotation_link where compare_id = $compare_id";
    my @results = &do_sql_2D($dbproc, $query);
    
    my %all_TUs;
    foreach my $result (@results) {
        my ($gene_id, $model_id) = @$result;
        $model_to_TU{$model_id} = $gene_id;
        $all_TUs{$gene_id} = 1;
    }
    
    my $total_number_genes_linked = scalar (keys %all_TUs);
    
    ## Determine how many genes were linked to FL assemblies:
    my $query = "select al.gene_id, al.cdna_acc "
        . " from annotation_link al, align_link al2, cdna_info ci "
        . " where al.compare_id = $compare_id and al.cdna_acc = al2.align_acc and al2.cdna_info_id = ci.id and ci.is_assembly = 1 and ci.is_fli = 1 ";
    my @results = &do_sql_2D($dbproc, $query);
    my %genes_linked_to_FL_assemblies;
    my %FL_assemblies_linked_to_genes;
    foreach my $result (@results) {
        my ($gene_id, $cdna_acc) = @$result;
        $genes_linked_to_FL_assemblies{$gene_id} = 1;
        $FL_assemblies_linked_to_genes{$cdna_acc} = 1;
    }
    my $number_genes_linked_to_FL_assemblies = scalar (keys %genes_linked_to_FL_assemblies);
    my $number_FL_assemblies_linked_to_genes = scalar (keys %FL_assemblies_linked_to_genes);
    
    my @genes_supported_by_PASA_assemblies;
    my @genes_supported_by_FL_assemblies;
    my @genes_supported_by_EST_assemblies_only;
    
    
    my @genes_failed_by_PASA_assemblies;
    my @genes_failed_by_FL_assemblies;
    my @genes_failed_by_EST_assemblies_only;
    
    
    ################# Analyze Successes ##############################
    
    my $root_query = "select distinct gene_id "
        . " from annotation_link al, status_link sl, status s, align_link al2, cdna_info ci "
        . " where al.cdna_acc = sl.cdna_acc and sl.status_id = s.status_id and sl.compare_id = $compare_id and al.compare_id = $compare_id and sl.cdna_acc = al2.align_acc and al2.cdna_info_id = ci.id and ci.is_assembly = 1 ";
    
    my $query = $root_query;
    $query .= "and (s.requires_update = 1 or (s.requires_update = 0 and s.fails_incorporation = 0)) ";
    my @results = &do_sql_2D($dbproc, $query);
    foreach my $result (@results) {
        my $gene = $result->[0];
        push (@genes_supported_by_PASA_assemblies , $gene);
    }
    
    my %fl;
    ## get those considered to be FL-only:
    $query .= " and ci.is_fli = 1";
    my @results = &do_sql_2D($dbproc, $query);
    foreach my $result (@results) {
        my $gene = $result->[0];
        $fl{$gene} = 1;
        push (@genes_supported_by_FL_assemblies, $gene);
    }
    
    ## find those genes that are not supported by FL assemblies:
    foreach my $gene (@genes_supported_by_PASA_assemblies) {
        unless ($fl{$gene}) {
            push (@genes_supported_by_EST_assemblies_only, $gene);
        }
    }
    
    my %genes_supported_by_PASA_assemblies_hash;
    my %genes_supported_by_FL_assemblies_hash;
    foreach my $gene (@genes_supported_by_PASA_assemblies) {
        $genes_supported_by_PASA_assemblies_hash{$gene} =1;
        if ($fl{$gene}) {
            $genes_supported_by_FL_assemblies_hash{$gene} = 1;
        }
    }
    
    ################## Analyze Failures #################################
    $query = $root_query;
    $query .= " and s.fails_incorporation = 1";
    
    my @genes_failed_but_successful_with_others;
    
    
    my @results = &do_sql_2D($dbproc, $query);
    foreach my $result (@results) {
        my $gene = $result->[0];
        push (@genes_failed_by_PASA_assemblies, $gene);
        if ($genes_supported_by_PASA_assemblies_hash{$gene}) {
            push (@genes_failed_but_successful_with_others, $gene);
        }
    }
    
    ## get those only FL
    my @genes_failed_but_successful_with_other_FLs;
    %fl = (); #reinit
    $query .= " and ci.is_fli = 1";
    my @results = &do_sql_2D($dbproc, $query);
    foreach my $result (@results) {
        my $gene = $result->[0];
        $fl{$gene} = 1;
        push (@genes_failed_by_FL_assemblies, $gene);
        if ($genes_supported_by_FL_assemblies_hash{$gene}) {
            push (@genes_failed_but_successful_with_other_FLs, $gene);
        }
    }
    
    ## get those failed w/ ESTs only:
    foreach my $gene (@genes_failed_by_PASA_assemblies) {
        unless ($fl{$gene}) {
            push (@genes_failed_by_EST_assemblies_only, $gene);
        }
    }
    
    ## get counts of assemblies supporting gene models:
    my $query = "select count(ci.cdna_acc) from "
        . " cdna_info ci, align_link al, status_link sl, status s "
        . " where ci.is_assembly = 1 and ci.is_fli = 1 and ci.id = al.cdna_info_id and al.align_acc = sl.cdna_acc and sl.compare_id = $compare_id and sl.status_id = s.status_id and (s.requires_update = 1 or (s.requires_update = 0 and s.fails_incorporation = 0))";
    my $number_FL_assemblies_supporting_models = &very_first_result_sql($dbproc, $query);
    
    
    my $query = "select count(ci.cdna_acc) "
        . " from cdna_info ci, align_link al, status_link sl, status s "
        . " where ci.is_assembly = 1 and ci.is_fli = 0 and ci.id = al.cdna_info_id and al.align_acc = sl.cdna_acc and sl.compare_id = $compare_id and sl.status_id = s.status_id and (s.requires_update = 1 or (s.requires_update = 0 and s.fails_incorporation = 0))";
    my $number_EST_assemblies_supporting_models = &very_first_result_sql($dbproc, $query);
    
    
    ##########################################
    #### Summarize DATA ######################
    ##########################################
    
    
    
    my $html_text =
        
        "<h2>Gene Comparison Overview</h2>\n"
        . "<table class='geneComparisonOverviewTable'>\n"
        . "<tr><th colspan=2>A total of <b>$total_number_genes_linked</b> genes mapped to PASA assemblies</th></tr>\n"
        . "<tr><td colspan=2><b>$number_FL_assemblies_linked_to_genes</b> FL-assemblies mapped to <b>$number_genes_linked_to_FL_assemblies</b> genes</td></tr>"
        . "<tr><th colspan=2 >Successful Incorporations by Genes</th></tr>\n"
        . "<tr><td>Genes supported by PASA assemblies</td><td><b>" . scalar (@genes_supported_by_PASA_assemblies) . "</b></td></tr>\n"
        . "<tr><td>Genes supported by <b>$number_FL_assemblies_supporting_models</b> FL assemblies</td><td><b>" . scalar (@genes_supported_by_FL_assemblies) . "</b></td></tr>\n"
        . "<tr><td>Genes supported by <b>$number_EST_assemblies_supporting_models</b> EST assemblies only <i>(no FL assemblies)</i></td><td><b>" . scalar (@genes_supported_by_EST_assemblies_only) . "</b></td></tr>\n"
        . "<tr><th colspan=2>Failed Incorporations by Genes</th></tr>\n"
        . "<tr><td>Genes failing PASA assembly incorporation</td><td><b>" . scalar (@genes_failed_by_PASA_assemblies) . "</b></td></tr>\n"
        . "<tr><td colspan=2><i>of these, <b>" . scalar (@genes_failed_but_successful_with_others) . "</b> genes successfully incorporate other PASA assemblies</i></td></tr>\n"
        . "<tr><td>Genes failing FL assembly incorporation</td><td><b>" . scalar (@genes_failed_by_FL_assemblies) . "</b></td></tr>\n"
        . "<tr><td colspan=2><i>but <b>" . scalar (@genes_failed_but_successful_with_other_FLs) . "</b> of these genes incorporate other FL assemblies successfully.</td></tr>\n"
        . "<tr><td>Genes failing EST assemblies only</td><td><b>" . scalar (@genes_failed_by_EST_assemblies_only) . "</b></td></tr>\n"
        . "</table>\n";
    $html_text .= "<i>The same gene may be counted in multiple categories.</i>\n";
    
    return ($html_text);
}


