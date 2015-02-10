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


unless ($db) {
    die "Must set the db parameter.\n";
}


my $css_common_text = &Pasa_CGI::get_CSS($0);    
print $cgi->start_html(-title => "Alternative Splicing Report for $db",
                       -head => style( { type => "text/css" }, $css_common_text ),   
    );

my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $mysql_ro_user = &Pasa_conf::getParam("MYSQL_RO_USER");
my $mysql_ro_password = &Pasa_conf::getParam("MYSQL_RO_PASSWORD");

my ($dbproc) = &connect_to_db($mysql_server,$db,$mysql_ro_user,$mysql_ro_password);

my %tokens;
my $query = "select token_id, alt_splice_token from alt_splice_tokens";
my @results = &do_sql_2D($dbproc, $query);
foreach my $result (@results) {
    my ($token_id, $alt_splice_token) = @$result;
    $tokens{$token_id} = $alt_splice_token;
}


if (my $token_id = $cgi->param('token_id') ) {
    &generate_assembly_listing( { 
        datatype => "token_id",
        token_id => $token_id
            }
                                );
    
}

elsif (my $alt_splice_type = $cgi->param('type')) {
    &generate_assembly_listing( { 
        datatype => "type",
        type => $alt_splice_type, 
        subtype => $cgi->param('subtype'), 
    }
                                );
}

else {
    &generate_report_stats();
}

print $cgi->end_html();
exit(0);


####
sub generate_report_stats() {
    

   
    
    my $max_compare_id = &Ath1_cdnas::get_max_compare_id($dbproc) || 0;
    
    print "<h1>Alternative Splicing Summary for $db</h1>\n";
    
    ##################################
    ## total number of transcript involved in alternative splicing:
    print "<h2>PASA assembly statistics</h2>\n"
        . "<i>*assemblies may be assigned to multiple categories</i>\n";
    print "<table border=1>\n";
    print "<tr><th>PASA assembly attribute</th><th># assemblies</th><th>% assemblies</th><th># subclusters</th><th>% subclusters</th><th># genes</th><th>% genes</th></tr>\n";
    
    
    my $query = "select count(*) from cdna_info where is_assembly = 1 ";
    my $total_number_assemblies = &very_first_result_sql($dbproc, $query);
    

    $query = "select count(*) from alt_splice_token_assignment";
    my $total_transcripts_alt_splicing = &very_first_result_sql($dbproc, $query);
    
    my $percentage_total_assemblies = "NA";
    if ($total_number_assemblies != 0) {
        $percentage_total_assemblies = sprintf ("%.2f", $total_transcripts_alt_splicing / $total_number_assemblies * 100);
    }
    
    ## number of genes involved in alternative splicing:
    
    $query = "select annotation_version from annotation_compare where compare_id = $max_compare_id";
    my $annot_version = &very_first_result_sql($dbproc, $query) || 0;
    
    $query = "select count(distinct gene_id) from annotation_store where annotation_version = $annot_version";
    my $total_number_genes = &very_first_result_sql($dbproc, $query);
    
    $query = "select count(distinct al.gene_id) from annotation_link al, alt_splice_token_assignment asta where al.cdna_acc = asta.cdna_acc and al.compare_id = $max_compare_id ";
    my $total_genes_alt_spliced = &very_first_result_sql($dbproc, $query);
    
    my $percentage_total_genes = "NA";
    if ($total_number_genes != 0) {
        $percentage_total_genes = sprintf ("%.2f", $total_genes_alt_spliced / $total_number_genes * 100);
    }

    ## number of subclusters alt spliced
    
    $query = "select count(distinct subcluster_id) from subcluster_link";
    my $total_number_subclusters = &very_first_result_sql($dbproc, $query);
    
    $query = "select count(distinct subcluster_id) from subcluster_link sl, alt_splice_token_assignment asta where sl.cdna_acc = asta.cdna_acc";
    my $number_subclusters_alt_spliced =  &very_first_result_sql($dbproc, $query);
    
    my $percentage_total_subclusters = "NA";
    if ($total_number_subclusters != 0) {
        $percentage_total_subclusters = sprintf ("%.2f", $number_subclusters_alt_spliced / $total_number_subclusters * 100);
    }

    
    print "<tr><td>Counts involved in alt-splicing:</td>"
        . "<td>$total_transcripts_alt_splicing</td><td>$percentage_total_assemblies</td>"
        . "<td>$number_subclusters_alt_spliced</td><td>$percentage_total_subclusters</td>"
        . "<td>$total_genes_alt_spliced</td><td>$percentage_total_genes</td></tr>\n";
    
    


    ## get alt-splice type:
    my @alt_splice_types;
    $query = "select distinct type from splice_variation";
    my @results = &do_sql_2D($dbproc, $query);
    foreach my $result (@results) {
        my ($type) = (@$result);
        push (@alt_splice_types, $type);
    }
    
    @alt_splice_types = sort (@alt_splice_types);
    foreach my $alt_splice_type (@alt_splice_types) {

        ## count assemblies:
        my $query = "select count(distinct cdna_acc) from splice_variation where type = ?";
        my $count_assemblies = &very_first_result_sql($dbproc, $query, $alt_splice_type);

        my $percentage_assemblies = "NA";
        if ($total_transcripts_alt_splicing != 0) {
            $percentage_assemblies = sprintf ("%.2f", $count_assemblies / $total_transcripts_alt_splicing * 100);
        }

        ## count subclusters:
        $query = "select count(distinct subcluster_id) from subcluster_link sl, splice_variation sv where sv.type = ? and sv.cdna_acc = sl.cdna_acc";
        my $count_subclusters = &very_first_result_sql($dbproc, $query, $alt_splice_type);
        
        my $percentage_subclusters = "NA";
        if ($number_subclusters_alt_spliced != 0) {
            $percentage_subclusters = sprintf ("%.2f", $count_subclusters / $number_subclusters_alt_spliced * 100);
        }


        ## count genes:
        $query = "select count(distinct al.gene_id) from annotation_link al, splice_variation sv where sv.type = ? and sv.cdna_acc = al.cdna_acc and al.compare_id = ?";
        my $count_genes = &very_first_result_sql($dbproc, $query, $alt_splice_type, $max_compare_id);
        
        my $percentage_genes = "NA";
        if ($total_genes_alt_spliced != 0) {
            $percentage_genes = sprintf ("%.2f", $count_genes / $total_genes_alt_spliced * 100);
        }
        
        print "<tr><td>$alt_splice_type</td><td><a href=\"alt_splice_report.cgi?db=$db&type=$alt_splice_type\">$count_assemblies</td><td>$percentage_assemblies</td><td>$count_subclusters</td><td>$percentage_subclusters</td><td>$count_genes</td><td>$percentage_genes</td></tr>\n";
    }
    print "</table>\n";
    
    print "<p id=\"altSpliceSearch\"><a href=\"search_alt_splice.cgi?db=$db\">search</a> for alt-splice variants</p>\n";   
 
    #########################  Individual feature count:
    #########################
    print "<h2>Counts of individual features</h2>\n";
    print "<table border=1>\n";
    print "<tr><th>feature</th><th>count</th></tr>\n";
    
    foreach my $type ("alt_acceptor", "alt_donor", "retained_intron", 
                      "spliced_intron", 
                      "retained_exon", 
                      "skipped_exon", 
                      "alternate_exon", "starts_in_intron", "ends_in_intron") {
        
        my $query = "select count(distinct c.annotdb_asmbl_id, s.lend, s.rend, s.orient) "
            . " from splice_variation s, align_link al, clusters c  "
            . " where type = ? and s.cdna_acc = al.align_acc and al.cluster_id = c.cluster_id ";
        my $count = &very_first_result_sql($dbproc, $query, $type);
        print "<tr><td>number of <b>$type</b>s</td><td>$count</td></tr>\n";
    }
    print "</table>\n";
    

    #################  Distribution of subfeatures within these events:
    
    foreach my $type ("skipped_exon", "alternate_exon") {
        
        my %subfeature_counts;
        my %seen_exon;
        my $query = "select c.annotdb_asmbl_id, s.lend, s.rend, s.orient, s.num_subfeatures_included "
            . " from splice_variation s, align_link al, clusters c "
            . " where type = \"$type\" and s.cdna_acc = al.align_acc and al.cluster_id = c.cluster_id "
            . " order by num_subfeatures_included desc";
        
        my @results = &do_sql_2D($dbproc, $query);
        foreach my $result (@results) {
            my ($asmbl_id, $lend, $rend, $orient, $num_subfeatures) = @$result;
            my $feature_key = "$asmbl_id,$lend,$rend,$orient";
            if ($seen_exon{$feature_key}) {
                next;
            }
            $seen_exon{$feature_key} = 1;
            $subfeature_counts{$num_subfeatures}++;
        }
        
        print "<h2>Distribution of exon counts within $type events</h2>\n"
            . "<table border=1>\n"
            . "<tr><th>num exons in $type</th><th>num $type events</th>\n";
        foreach my $num_subfeats (reverse sort {$a<=>$b} keys %subfeature_counts) {
            my $num_subfeat_count = $subfeature_counts{$num_subfeats};
            print "<tr><td>$num_subfeats</td>"
                . "<td><a href=\"search_alt_splice.cgi?db=$db&SEARCH=1&spliceVariationTypes=$type&numSubVariations=$num_subfeats\">$num_subfeat_count</a></td></tr>\n";
        }
        print "</table>\n";
    }

    my $query = "select count(distinct c.annotdb_asmbl_id, s.lend, s.rend, s.orient) "
        . " from splice_variation s, align_link al, clusters c "
        . " where type = 'retained_exon' and subtype = 'alternate_internal_exons' "
        . " and s.cdna_acc = al.align_acc and al.cluster_id = c.cluster_id ";
    my $count_alternate_internal_exons = &very_first_result_sql($dbproc, $query);

    print "<p><a href=\"alt_splice_report.cgi?db=$db&type=retained_exon&subtype=alternate_internal_exons\">$count_alternate_internal_exons</a> internal exons appear to be alternate internal exons.\n";
    
    
    ## Links to graphs
    print "<h2>Graphs</h2>\n"
        . "<ul>"
        . "<li><a href=\"variation_support.cgi?db=$db&gradation=0.1\" target=splice_variation_support_$db >Relative support for individual variations</a>\n"
        . "<li><a href=\"alt_donor_acceptor_deltas.cgi?db=$db&gradation=1&max_delta=100\" target=_delta_splice_diffs_$db\">Lengths between alternate splice sites.</a>"
        . "</ul>\n";
    
    ## Links to additional resources
    print "<h2>Additional analyses</h2>\n"
        . "<ul>\n"
        . "<li><a href=\"FL_ref_assembly_compare.cgi?db=$db\">localizations to UTRs or CDS regions using a FL-assembly based reference gene structure.</a>\n"
        . "</ul>\n";
    
    
    


    #############################
    ## Combinations
    
    print "<h2>Combinations of alternative splicing variations</h2>\n";
    print "<i>*categories are mutually exclusive</i>\n";
    print "<table border=1>\n"
        . "<tr><th>alt-splice combination</th><th>count</th></tr>\n";
    
    $query = "select token_id, count(*) from alt_splice_token_assignment group by token_id order by 2 desc";
    @results = &do_sql_2D($dbproc, $query);
    foreach my $result (@results) {
        my ($token_id, $count) = @$result;
        print "<tr><td>$tokens{$token_id}</td><td><a href=\"alt_splice_report.cgi?db=$db&token_id=$token_id\">$count</td></tr>\n";
    }
    print "</table>\n";
    
}

    
####
sub generate_assembly_listing {
    my ($data_href) = @_; #$data_type, $value) = @_;
    
    my $data_type = $data_href->{datatype};
    my $value = $data_href->{$data_type};
    

    my $query = "";
    my @params;

    if ($data_type eq "token_id") {
        $query = "select cdna_acc from alt_splice_token_assignment where token_id = ?";
        @params = ($value);
        print "<h2>Assemblies assigned to $tokens{$value}</h2>\n";
    }
    
    elsif ($data_type eq "type") {
    
        my $subtype = $data_href->{subtype};
    
        $query = "select distinct cdna_acc from splice_variation where type = ?";
        @params = ($value);
        
        if ($subtype) {
            $query .= " and subtype = ? ";
            push (@params, $subtype);

            $value .= ", subtype: $subtype";
        }
        
        print "<h2>Alternative splicing with classification of [$value]</h2>\n";
    }
    else {
        die "Error, no data type specified";
    }
    
    print "<table border=1>\n"
        . "<tr><th>count</th><th>pasa assembly accession</th></tr>\n";
    
    my @results = &do_sql_2D($dbproc, $query, @params);
    
    my $count = 0;
    foreach my $result(@results) {
        $count++;
        my ($cdna_acc) = @$result;
        print "<tr><td>$count</td><td><a href=\"assembly_alt_splice_info.cgi?db=$db&cdna_acc=$cdna_acc\">$cdna_acc</td></tr>\n";
    }
    print "</table>\n";
    
    
}



    
    
    
    
