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
use POSIX qw(ceil);

$|++;
my $cgi = new CGI;
print $cgi->header('text/html');
my $db = $cgi->param('db');

my $show_entries_flag = $cgi->param('SHOW_ENTRIES');
my $bin_size = $cgi->param('gradation') || 1;
my $max_delta = $cgi->param('max_delta') || 0;


our $SEE = $cgi->param('SEE') || $cgi->param('DEBUG');

print $cgi->start_html(-title=> "Lengths Between Alternative Acceptors and Donors" );

unless ($ENV{WEBSERVER_TMP}) {
    $ENV{WEBSERVER_TMP} = "/tmp";
}

unless ($db) {
    die "Must set the db parameter.\n";
}


my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $mysql_ro_user = &Pasa_conf::getParam("MYSQL_RO_USER");
my $mysql_ro_password = &Pasa_conf::getParam("MYSQL_RO_PASSWORD");

my ($dbproc) = &connect_to_db($mysql_server,$db,$mysql_ro_user,$mysql_ro_password);

my $url = "alt_donor_acceptor_deltas.cgi?db=$db&bin_size=$bin_size&max_delta=$max_delta";
if ($show_entries_flag) {
    print "<a href=\"$url&SHOW_ENTRIES=0\">[hide text data]</a>\n";
}
else {
    print "<a href=\"$url&SHOW_ENTRIES=1\">[show text data]</a>\n";
}

## allow for max deltas:

$url = "alt_donor_acceptor_deltas.cgi?db=$db&bin_size=$bin_size&SHOW_ENTRIES=$show_entries_flag";
print " [maximum delta = ";
foreach my $i (10, 20, 50, 100, 200, "none") {
    
    my $val = $i;
    if ($val eq "none") {
        $val = "";
    } else {
        $val = "&max_delta=$i";
    }
    
    print "<a href=\"$url&$val\">$i</a>";
    if ($i ne "none") {
        print ", ";
    }
}
print "]\n";


print "<hr>\n";

my %genome_contig;

if ($show_entries_flag) {
    ## get genomic contig info
    my $query = "select c.annotdb_asmbl_id, al.align_acc from clusters c, align_link al, cdna_info ci where c.cluster_id = al.cluster_id and al.cdna_info_id = ci.id and ci.is_assembly = 1";
    my @results = &do_sql_2D($dbproc, $query);
    foreach my $result (@results) {
        my ($contig_id, $cdna_acc) = @$result;
        $genome_contig{$cdna_acc} = $contig_id;
    }
}




foreach my $splice_type ("alt_acceptor", "alt_donor") {


    # put named anchors to help navigate the data
    if ($show_entries_flag) { 
        if ($splice_type eq "alt_donor") {
            print "<a name='alt_donor'>&nbsp;</a>\n";
            print "<a href='#alt_acceptor'>[goto alt-acceptor data]</a>\n";
        } else {
            print "<a name='alt_acceptor'>&nbsp;</a>\n";
            print " <a href=\"#alt_donor\">[goto alt-donor data]</a>\n";
        }
    }
    
    my $query = qq {select sv1.cdna_acc, sv1.orient, sv1.lend, sv1.rend, sv2.cdna_acc, sv2.lend, sv2.rend 
                        from splice_variation sv1, splice_variation sv2, alt_splice_link asl
                        where sv1.sv_id = asl.sv_id_A 
                        and asl.sv_id_B = sv2.sv_id
                        and asl.sv_id_A < asl.sv_id_B 
                        and sv1.type = "$splice_type"
                        and sv2.type = "$splice_type"
                        
                    };
    
    

    my @structs;
    
    my %counts;

    my @results = &do_sql_2D($dbproc, $query);
    
    my $data_text = "<pre><b>splice_type\tcontig\torientation\taccA\tlendA\trendA\taccB\tlendB\trendB\tdelta</b>\n";
    
    foreach my $result (@results) {
        my ($cdna_acc_A, $orient, $lend_A, $rend_A, $cdna_acc_B, $lend_B, $rend_B) = @$result;
        
        my $delta = abs ($lend_A - $lend_B);
        
        my $bin_no = ceil($delta/$bin_size);
        $counts{$bin_no}++;
        
        if ($show_entries_flag && ($max_delta == 0 || $delta <= $max_delta)) {
            my $contig = $genome_contig{$cdna_acc_A};
            $data_text .= "$splice_type\t$contig\t$orient\t$cdna_acc_A\t$lend_A\t$rend_A\t$cdna_acc_B\t$lend_B\t$rend_B\t$delta\n";
        }
        
    }
    $data_text .= "</pre>\n";
    


    my @data = ([], []);
    foreach my $bin (sort {$a<=>$b} keys %counts) {
        my $x_value = $bin * $bin_size;
        my $y_value = $counts{$bin};
        
        if ($max_delta==0 || $x_value <= $max_delta) {
            push (@{$data[0]}, $x_value);
            push (@{$data[1]}, $y_value);
        }
    }
    
    
    my $image_file = "$ENV{WEBSERVER_TMP}/$splice_type.$$.png";
    my $graph = new GD::Graph::lines (800,300);
    $graph->set ( x_label => "delta in basepairs",
                  y_label => "number of alternate splice sites",
                  title => "deltas for $splice_type",
                  x_labels_vertical => 1,
                  show_values => 1,
                  #bar_spacing => 4,
                  dclrs => ["blue"],
                  #x_max_value => 1,
                  );
    
    my $gd = $graph->plot (\@data);
    
    open (my $fh, ">$image_file") or die $!;
    binmode $fh;
    print $fh $gd->png;
    close $fh;
    
    print "<p><img src=\"show_png.cgi?image=$image_file\">";
    
    if ($show_entries_flag) {
        print $data_text;
    }
    
    print "<hr>\n";
    
}


print $cgi->end_html();

exit(0);


