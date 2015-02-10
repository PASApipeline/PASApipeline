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
use GD::Graph::lines;
use GD::Graph::bars;
use POSIX qw(ceil);

$|++;
my $cgi = new CGI;
print $cgi->header('text/html');
my $db = $cgi->param('db');

my $show_entries_flag = $cgi->param('SHOW_ENTRIES');
my $gradation = $cgi->param('gradation') || 0.1;
my $exclude_1_1 = $cgi->param('exclude_1_1') || 0;


our $SEE = $cgi->param('SEE') || $cgi->param('DEBUG');

print $cgi->start_html(-title=> "splicing variation support" );

unless ($ENV{WEBSERVER_TMP}) {
    $ENV{WEBSERVER_TMP} = "/tmp";
}

unless ($db) {
    die "Must set the db parameter.\n";
}



print "<a name='top'>&nbsp;</a>";
print "<p>PASA assemblies demonstrating splicing variations were examined to identify the underlying ESTs and cDNAs responsible for the variation.  The ratio between the number of underlying transcripts demonstrating the splicing variation is used as an indicator of the relative support for the variation.\n<hr>\n";


my @variations = qw (
                     alt_acceptor
                     alt_donor
                     retained_intron
                     alternate_exon
                     skipped_exon
                     starts_in_intron
                     ends_in_intron
                     );


## must divide symmetric variations by 2 if storing a-b and b-a links.
my %symmetric_variation = ( alt_donor => 2,
                            alt_acceptor => 2,
                            alternate_exon => 2,
                            );

my $url = "variation_support.cgi?db=$db&gradation=$gradation&exclude_1_1=$exclude_1_1";
if ($show_entries_flag) {
    print "<a href=\"$url&SHOW_ENTRIES=0\">[hide text data]</a><br>\n";
        
} else {
    print "<a href=\"$url&SHOW_ENTRIES=1\">[show text data]</a><br>\n";
}
 
$url = "variation_support.cgi?db=$db&gradation=$gradation&SHOW_ENTRIES=$show_entries_flag";
if ($exclude_1_1) {
    print " [<a href=\"$url&exclude_1_1=0\">Include</a> cases where only 1 transcript supports each variation]<br>\n";
} else {
    print " [<a href=\"$url&exclude_1_1=1\">Exclude</a> cases where only 1 transcript supports each variation]<br>\n";
}


my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $mysql_ro_user = &Pasa_conf::getParam("MYSQL_RO_USER");
my $mysql_ro_password = &Pasa_conf::getParam("MYSQL_RO_PASSWORD");

my ($dbproc) = &connect_to_db($mysql_server,$db,$mysql_ro_user,$mysql_ro_password);

my %genome_contig;

if ($show_entries_flag) {
    ## get genomic contig info
    my $query = "select c.annotdb_asmbl_id, cl.cdna_acc from clusters c, cluster_link cl where c.cluster_id = cl.cluster_id and cl.is_assembly = 1";
    my @results = &do_sql_2D($dbproc, $query);
    foreach my $result (@results) {
        my ($contig_id, $cdna_acc) = @$result;
        $genome_contig{$cdna_acc} = $contig_id;
    }
}


foreach my $variation (@variations) {
    
    ## include named anchors
    print "goto: ";
    foreach my $variation ("top", @variations) {
        print "<a href=\"\#$variation\">$variation,</a>\n";
    }
    print "<a name=\"$variation\">&nbsp;</a>\n";
    print "<br>\n";
    
    my $query = qq { select sv.cdna_acc, svs.cdna_acc, svs.num_transcripts_A, svs.num_transcripts_B 
                         from splice_variation sv, splice_variation_support svs
                         where sv.sv_id = svs.sv_id 
                         and sv.type = ? };
    my @results = &do_sql_2D($dbproc, $query, $variation);

        
    my @entries;
    my $max_ratio = 0;
    foreach my $result (@results) {
        my ($cdna_acc, $other_acc, $num_transcripts_A, $num_transcripts_B) = @$result;
        
        my ($x, $y) = ($num_transcripts_A, $num_transcripts_B);

        if ($x == 0 || $y == 0) { next; } #don't count cases where no transcripts found as support
        

        if ($exclude_1_1 && $x == 1 && $y == 1) {
            next; #excluding entry
        }
        

        if ($symmetric_variation{$variation}) {
            ($x, $y) = sort {$a<=>$b} ($num_transcripts_A, $num_transcripts_B);
        }
        my $ratio = $x / $y;

        push (@entries, { ratio => $ratio,
                          cdna_acc => $cdna_acc,
                          num_A => $num_transcripts_A,
                          num_B => $num_transcripts_B,
                          other_acc => $other_acc,
                      } );

    
        if ($ratio > $max_ratio) {
            $max_ratio = $ratio;
        }
    }

    @entries = reverse sort { $a->{ratio} <=> $b->{ratio} } @entries;

    ## create graph
    
    ## bin at every gradation
    my %counts;
    
    my $graph;
    my @data = ([], []);
    my $totals = 0;
    
    if (my $divider = $symmetric_variation{$variation}) {
        # init counts:
        for (my $i = $gradation; $i <= ceil($max_ratio); $i += $gradation) {
            $counts{$i} = 0;
        }
        $graph = new GD::Graph::bars (500,300);
        
        ## map to integer values
        my $grad_bin = int($gradation * 100 + 0.5); 
        foreach my $entry (@entries) {
            my $ratio = $entry->{ratio};
            my $ratio_adj = int ($ratio * 100 + 0.5); # round up
            
            my $bin = int (($ratio_adj -1) / $grad_bin);
            $bin++;
            if ($bin == 0) {
                $bin = 1;
            }
            
            my $bin_val = $bin * $gradation;
            
            $counts{$bin_val}++;
        }
        
        foreach my $datapoint (sort {$a<=>$b} keys %counts) {
            my $value = $counts{$datapoint};
            
            $value /= $divider; # adjust for symmetry (fully symmetric entries are included twice)
            $value = int($value);
            
            push (@{$data[0]}, $datapoint); # x-value
            
            push (@{$data[1]}, $value); # y-value
            $totals += $value;
        }
        
    
        
        
        $graph->set ( x_label => "ratio of transcript support",
                      y_label => "number of assembly comparisons",
                      title => "varied support for type [$variation]",
                      x_labels_vertical => 1,
                      show_values => 1,
                      bar_spacing => 4,
                      dclrs => ["green"],
                      #x_all_ticks => 1,
                      #x_tick_number => "auto",
                      );
        

    }else {
        ## symmetric variation
        $graph = new GD::Graph::bars (500,300);
        my $counts_below = 0;
        my $counts_at_1 = 0;
        my $counts_above = 1;
        foreach my $entry (@entries) {
            my $ratio = $entry->{ratio};
            if ($ratio < 1) {
                $counts_below++;
            }
            elsif ($ratio > 1) {
                $counts_above++;
            }
            else {
                $counts_at_1++;
            }
            $totals++;
        }
        push (@{$data[0]}, "<", "1", ">");
        push (@{$data[1]}, $counts_below, $counts_at_1, $counts_above);


        $graph->set ( x_label => "ratio of transcript support",
                      y_label => "number of assembly comparisons",
                      title => "varied support for type [$variation]",
                      #x_labels_vertical => 1,
                      show_values => 1,
                      bar_spacing => 4,
                      dclrs => ["red"],
                      #x_all_ticks => 1,
                      #x_tick_number => "auto",
                      );


    }
    
    my $gd = $graph->plot (\@data);
    my $image_file = "$ENV{WEBSERVER_TMP}/$variation.$$.png";
    open (my $fh, ">$image_file") or die $!;
    binmode $fh;
    print $fh $gd->png;
    close $fh;

    print "<h2>[$variation] Splicing Variation Support</h2>\n";
    print "<p><img src=\"show_png.cgi?image=$image_file\">";
    print "<br>total: $totals pairwise comparisons between pasa assemblies.\n";

    if ($show_entries_flag) {
        ## show data:
        print "<pre><b>contig\tcdna_accA\tcdna_accB\t#subset_transcipts_A\t#subset_transcripts_B\tratio_support</b>\n";
        foreach my $entry (@entries) {
            my ($cdna_acc, $other_acc, $num_A, $num_B, $ratio) = ($entry->{cdna_acc},
                                                                  $entry->{other_acc},
                                                                  $entry->{num_A},
                                                                  $entry->{num_B},
                                                                  $entry->{ratio},
                                                                  );
            my $contig = $genome_contig{$cdna_acc};
            print "$contig\t$cdna_acc\t$other_acc\t$num_A\t$num_B\t$ratio\n";
        }
        print "</pre>\n";
    }
    
    print "<hr>\n";
    
}


print $cgi->end_html();

exit(0);


