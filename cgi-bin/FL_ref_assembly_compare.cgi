#!/usr/bin/env perl

use strict;
use warnings;
use Pasa_init;
use Pasa_conf;
use DBI;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Data::Dumper;
use Mysql_connect;
use Ath1_cdnas;
use GD::Graph::lines;
use CGI::Pretty ":standard";
use Pasa_CGI;


$|++;
my $cgi = new CGI;
print $cgi->header('text/html');
my $db = $cgi->param('db');
our $SEE = $cgi->param('SEE') || $cgi->param('DEBUG');

my $css_common_text = &Pasa_CGI::get_common_CSS();


print $cgi->start_html(-title => "Using a FL-assembly reference gene structure to localize splicing variations ($db)",
                       -head => style( { type => "text/css" }, $css_common_text ),
                       );

unless ($db) {
    die "Must set the db parameter.\n";
}

unless ($ENV{WEBSERVER_TMP}) {
    $ENV{WEBSERVER_TMP} = "/tmp";
}

my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $mysql_ro_user = &Pasa_conf::getParam("MYSQL_RO_USER");
my $mysql_ro_password = &Pasa_conf::getParam("MYSQL_RO_PASSWORD");

my ($dbproc) = &connect_to_db($mysql_server,$db,$mysql_ro_user,$mysql_ro_password);


print "<h2>Localizing Splicing Variations to UTRs and CDS regions</h2>\n";
print "<p>The FL-assembly encoding the longest protein is chosen as a reference to which all other assemblies within each corresponding subcluster are compared</p>\n"; 


## get count of reference FL-cdnas

my $query = "select count(distinct template_cdna_acc) from alt_splice_FL_compare";
my $count_of_reference_FL_assemblies = &very_first_result_sql($dbproc, $query);

$query = "select count(*) from alt_splice_FL_compare";
my $count_of_variations = &very_first_result_sql($dbproc, $query);

print "<b>$count_of_variations</b> splicing variations were mapped to <b>$count_of_reference_FL_assemblies</b> reference FL-assemblies.\n";

## generate global statistics:

my %localization_counts;
my %localizations;

$query = "select type, classification from alt_splice_FL_compare";
my @results = &do_sql_2D($dbproc, $query);

foreach my $result (@results) {
    my ($type, $classification) = @$result;
    my @indiv_classes = split (/,/, $classification);
    foreach my $indiv_class (@indiv_classes) {
        $localization_counts{$type}->{$indiv_class}++;
        $localizations{$indiv_class}=1;
    }
}

## Statistics when comparing all subclustered assemblies to the reference FL-entry
print "<table border=1>\n";
foreach my $type (sort keys %localization_counts) {
    my $data_ref = $localization_counts{$type};
    my $total = 0;
    foreach my $localization (keys %$data_ref) {
        my $count = $data_ref->{$localization};
        $total += $count;
    }
    
    
    print "<tr><td>$type<br>$total events</td><td>"
        . "<table>\n";
   

    foreach my $localization (sort keys %localizations) {
        my $count = $data_ref->{$localization} || 0;
        my $percentage = sprintf ("%.2f", $count / $total * 100);
        print "<tr><td>$localization</td><td>$count</td><td>($percentage%)</td></tr>\n";
    }
    print "</table></td></tr>\n";
}
print "</table>\n";


## Statistics when comparing only FL-assemblies to the reference
print "<hr>\n<h2>Analysis of Splicing Variations Only Between FL-assemblies</h2>\n";
print "<p>Within a subcluster of FL-assemblies, comparisons are performed against the chosen reference FL-assembly</p>\n";

$query = "select count(distinct template_acc) from alt_splice_FL_to_FL_compare where same_frame_exists = 1 and diff_in_cds = 1";
my $number_ref_FL = &very_first_result_sql($dbproc, $query);

$query = "select count(other_acc) from alt_splice_FL_to_FL_compare where same_frame_exists = 1 and diff_in_cds = 1";
my $number_other_FL = &very_first_result_sql($dbproc, $query);

$query = "select sum(num_variations) from alt_splice_FL_to_FL_compare where same_frame_exists = 1 and diff_in_cds = 1";
my $total_variations = &very_first_result_sql($dbproc, $query);

$query = "select count(*) from alt_splice_FL_to_FL_compare where same_frame_exists = 1 and diff_in_cds = 1 and num_variations = 1";
my $total_variations_where_one_diff_only = &very_first_result_sql($dbproc, $query);

$query = "select count(*) from alt_splice_FL_to_FL_compare where same_frame_exists = 1 and diff_in_cds = 1 and num_variations > 1";
my $total_variations_where_multiple_diffs = &very_first_result_sql($dbproc, $query);




## Analyze Changes in Proteins resulting from variation by comparing only FL-entries
print "<p><b>$number_other_FL</b> FL entries were compared to <b>$number_ref_FL</b> reference FL entries, examining <b>$total_variations</b> splicing variations.  Of these, <b>$total_variations_where_one_diff_only</b> comparisons involved only one splicing variation, and <b>$total_variations_where_multiple_diffs</b> comparisons involved multiple variations.<br><i>Only comparisons between FL entries encoding protein with ORFs overlapping in the same frame were examined, and only cases where a difference in the cds region was found. The reference structure is required to encode a protein at least 200 residues long.  Also, since we want to analyze the impact of a given splice variation as it relates to the reference protein, the reference structure is devoid of the asymetric variations (ends_in_intron, starts_in_intron) as inferred from a comparison to the other FL-assembly.</i></p>";


# first, study cases where any variation exists.
$query = "select template_acc, other_acc, frame_change, percent_prot_length from alt_splice_FL_to_FL_compare where diff_in_cds = 1 and same_frame_exists = 1";
my @entries = &get_entries($query);
print "<h2>Changes in proteins due to splicing differences</h2><i> (multiple splicing variations may exist in single comparisons)</i><br>\n";
&generate_fl_fl_compare_report(\@entries, "all");


# now, study cases where only a single variation is allowed.
print "<hr>\n";

$query = "select template_acc, other_acc, frame_change, percent_prot_length from alt_splice_FL_to_FL_compare where diff_in_cds = 1 and same_frame_exists = 1 and num_variations = 1";
@entries = &get_entries($query);
print "<h2>Changes in proteins due to splicing differences, restricted to single variations</h2>\n";
&generate_fl_fl_compare_report(\@entries, "single");


print $cgi->end_html();
$dbproc->disconnect;

exit(0);


sub get_entries {

    my $query = shift;
    
    my @entries;
    @results = &do_sql_2D($dbproc, $query);
    foreach my $result (@results) {
        my ($template_acc, $other_acc, $frame_change, $percent_prot_length) = @$result;
        
        my $query = "select type from splice_variation sv, splice_variation_support svs where sv.sv_id = svs.sv_id and svs.cdna_acc = ? and sv.cdna_acc = ?";
        my @results = &do_sql_2D($dbproc, $query, $template_acc, $other_acc);
        
        my @types;
        foreach my $result (@results) {
            my ($type) = @$result;
            push (@types, $type);
        }
        
        my $struct = { template_acc => $template_acc,
                       other_acc => $other_acc,
                       frame_change => $frame_change,
                       percent_prot_length => $percent_prot_length,
                       variation_types => \@types,
                       
                   };
        
        push (@entries, $struct);
        
    }

    return (@entries);

}



####
sub generate_fl_fl_compare_report {
    
    my $entries_list_aref = shift;
    my $report_type = shift;

    my @entries = @$entries_list_aref;
    
    ## want:
    # counts of variation types
    # counts of how often a variation type changes the frame
    # average percentage of protein length

    my %data;
    # indexed by splicing type
    #  structured like so:
    #    $data{type}->{counts_total}
    #               ->{counts_changes_frame}
    #               ->{sum_percent_prot_length}
    
    
    my %type_to_bin_counter;

    foreach my $entry (@entries) {
        my $frame_change = $entry->{frame_change};
        my $percent_prot_length = $entry->{percent_prot_length};
        my $variation_types_aref = $entry->{variation_types};
        
        foreach my $type (@$variation_types_aref) {
            $data{$type}->{counts_total}++;
            $data{$type}->{sum_percent_prot_length} += $percent_prot_length;
            if ($frame_change) {
                $data{$type}->{counts_changes_frame}++;
            }

            &bin_type(\%type_to_bin_counter, $type, $percent_prot_length);
            
        }
    }
    
    ## generate summary
    
    print "<table border=1>\n";
    print "<tr><th>variation</th><th>total counts</th><th>counts frame changed</th><th>percentage change frame</th><th>average percentage of protein length</th></tr>\n";
    foreach my $type (sort keys %data) {
        my $stats_ref = $data{$type};
        print "<tr><td><a href=\"alt_splice_FL_FL_compare.cgi?db=$db&variation_type=$type&report=$report_type\">$type</a></td>"
            . "<td>" . ($stats_ref->{counts_total} || 0) . "</td>"
            . "<td>" . ($stats_ref->{counts_changes_frame} || 0) . "</td>"
            . "<td>" . sprintf ("%.2f", $stats_ref->{counts_changes_frame} / $stats_ref->{counts_total} * 100) . "</td>"
            . "<td>" . sprintf ("%.2f", $stats_ref->{sum_percent_prot_length} /  $stats_ref->{counts_total}) . "</td></tr>\n";
    }
    print "</table>\n";

    &draw_percent_prot_length_illustration(\%type_to_bin_counter, $report_type);
    
}


####
sub bin_type {
    my ($type_to_bin_counter_href, $type, $percent_prot_length) = @_;
    
    my $ratio = (int ($percent_prot_length / 20) ) * 20;
    
    $type_to_bin_counter_href->{$type}->{$ratio}++;
}


####
sub draw_percent_prot_length_illustration {
    my ($type_to_bin_counter_href, $report_type) = @_;

    my $img = "$ENV{WEBSERVER_TMP}/$report_type.prot_length.png";

    open (my $fh, ">$img") or die $!;
    binmode $fh;

    my $graph = GD::Graph::lines->new(600, 400);
    
    my @types;
    my @data;
    
    my %ratio_values;
    my %type_to_totals;

    foreach my $type (keys %$type_to_bin_counter_href) {
        push (@types, $type);
        my $data_href = $type_to_bin_counter_href->{$type};
        foreach my $ratio (keys %$data_href) {
            $ratio_values{$ratio} = 1;
            $type_to_totals{$type}+= $data_href->{$ratio}; ## count entries in this category
        }
    }
    
    my @ratios = sort {$a<=>$b} keys %ratio_values;
    
    #print "<pre>\n";
    #use Data::Dumper;
    #print Dumper ($type_to_bin_counter_href);


    my $max_percentage = 0;
    
    foreach my $type (@types) {
        my @data_points = ();
        
        my $total_for_type = $type_to_totals{$type};

        
        foreach my $ratio (@ratios) {
            
            my $number_of_type_at_ratio = $type_to_bin_counter_href->{$type}->{$ratio};
            if (defined $number_of_type_at_ratio) {
                
                my $normalized_percentage = sprintf ("%.1f", $number_of_type_at_ratio / $total_for_type * 100);
                #print "$type, x=$ratio, y=$normalized_percentage\n";
                push (@data_points, $normalized_percentage);
                if ($normalized_percentage > $max_percentage) {
                    $max_percentage = $normalized_percentage;
                }
            }
            else {
                push (@data_points, undef);
            }
        }
        push (@data, [@data_points]);
    }

    $max_percentage = int ($max_percentage + 0.5);
    $graph->set_legend(@types);
    $graph->set( 
                 x_label => '% of longest protein length',
                 y_label => '% of comparisons involving splicing class',
                 title => 'Protein Truncation Resulting from Splicing Variation',
                 y_max_value => $max_percentage + 5,
                 # skip_undef => 1,
                 transparent => 0,
                 );

    print $fh $graph->plot( [ \@ratios, @data] )->png(); 
    
    close $fh;

    print "<img src=\"show_png.cgi?image=$img\">";

}
        

    
