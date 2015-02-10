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
use Gene_obj;
use Gene_cdna_image;
use GD;
use Sequence_feature_drawer;
use TierFeatures;
use Storable qw (thaw);
use CGI::Pretty ":standard";
use Pasa_CGI;


my $cgi = new CGI();

print $cgi->header('text/html');

my $db = $cgi->param('db');
my $reference_cdna_acc = $cgi->param('reference_cdna_acc');
my $variation_type = $cgi->param('variation_type');
my $report = $cgi->param('report'); # all|single

unless ($db && $variation_type) {
    die "Need db and variation_type\n";
}

my $pid = $$;
unless ($ENV{WEBSERVER_TMP}) {
    $ENV{WEBSERVER_TMP} = "/tmp";
}

my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $mysql_ro_user = &Pasa_conf::getParam("MYSQL_RO_USER");
my $mysql_ro_password = &Pasa_conf::getParam("MYSQL_RO_PASSWORD");

my ($dbproc) = &connect_to_db($mysql_server,$db,$mysql_ro_user,$mysql_ro_password);

my $css_common_text = &Pasa_CGI::get_common_CSS();
print $cgi->start_html( -title => "Alt-splicing FL-FL compare $reference_cdna_acc ($db)",
                        -head => style( { type => "text/css" }, $css_common_text ),
                        );

if ($reference_cdna_acc) {
    &generate_cdna_report();
}
else {
    &generate_listing();
}

print $cgi->end_html();

$dbproc->disconnect;

exit(0);


####
sub generate_listing {
    
    print "<h2>FL comparisons to FL [$variation_type]</h2>\n";
    
    my $query = "select distinct template_acc, other_acc, frame_change, percent_prot_length, num_variations from alt_splice_FL_to_FL_compare asf, splice_variation sv, splice_variation_support svs where other_acc = sv.cdna_acc and type = \"$variation_type\" and diff_in_cds = 1 and same_frame_exists = 1 and svs.sv_id = sv.sv_id and svs.cdna_acc = asf.template_acc";
    
    if ($report eq "single") {
        $query .= " and num_variations = 1";
    }
    
    my @results = &do_sql_2D($dbproc, $query);
    
    print "<table border=1>\n"
        . "<tr><th>count</th><th>template_acc</th><th>other_acc</th><th>frame_change</th><th>percent protein length</th><th>number of variations</th></tr>\n";
    
    my $counter = 0;
    foreach my $result (@results) {
        my ($template_acc, $other_acc, $frame_change, $percent_prot_length, $num_variations) = @$result;
        
        $frame_change = ($frame_change) ? "Yes" : "No";
        
        $counter++;
        print "<tr><td>$counter</td><td><a href=\"alt_splice_FL_FL_compare.cgi?db=$db&variation_type=$variation_type&reference_cdna_acc=$template_acc\">$template_acc</a></td>"
            . "<td>$other_acc</td>"
            . "<td>$frame_change</td>"
            . "<td>$percent_prot_length</td>"
            . "<td>$num_variations</td></tr>\n";
    }
    print "</table>\n";
}


####
sub generate_cdna_report {
    my $query = "select asf.other_acc, asf.frame_change, asf.percent_prot_length, asf.num_variations, sv.type, sv.lend, sv.rend, sv.orient from splice_variation sv, alt_splice_FL_to_FL_compare asf, splice_variation_support svs where asf.other_acc = sv.cdna_acc and svs.sv_id = sv.sv_id and svs.cdna_acc = asf.template_acc and diff_in_cds = 1 and same_frame_exists = 1 and asf.template_acc = ?";
    
    my %other_acc_to_struct; #group by other_acc
    
    my @results = &do_sql_2D($dbproc, $query, $reference_cdna_acc);
    foreach my $result (@results) {
        my ($other_acc, $frame_change, $percent_prot_length, $num_variations, $type, $lend, $rend, $orient) = @$result;
        
        my $struct = $other_acc_to_struct{$other_acc};
        if (ref $struct) {
            ## add variation info to existing entry
            push (@{$struct->{variations}}, { type => $type,
                                              lend => $lend,
                                              rend => $rend,
                                              orient => $orient,
                                          } );
            
            
        } else {
            ## instantiate it
            $other_acc_to_struct{$other_acc} = { other_acc => $other_acc,
                                                 frame_change => $frame_change,
                                                 percent_prot_length => $percent_prot_length,
                                                 num_variations => $num_variations,
                                                 variations => [ 
                                                                 { type => $type,
                                                                   lend => $lend,
                                                                   rend => $rend,
                                                                   orient => $orient,
                                                               }
                                                                 ],
                                                 
                                             };
        }
        
    }
    
    
    print "<h2>FL-to-FL variation report for $reference_cdna_acc</h2>\n";

    foreach my $other_acc (keys %other_acc_to_struct) {
        my $struct = $other_acc_to_struct{$other_acc};
        
        my $variant_list_aref = $struct->{variations};
        
        print "<hr><h2>Comparison of FL-$reference_cdna_acc to FL-$other_acc</h2>\n";
        
        print "<table>\n"
            . "<tr><td>other accession</td><td><a href=\"assembly_alt_splice_info.cgi?db=$db&cdna_acc=$other_acc\">$other_acc</a></td></tr>\n"
            . "<tr><td>frame change?</td><td>" . (($struct->{frame_change}) ? "Yes" : "No") . "</td></tr>\n"
            . "<tr><td>percentage of protein length</td><td>" . $struct->{percent_prot_length} . "</td></tr>\n"
            . "<tr><td>number of variations</td><td>" . $struct->{num_variations} . "</td></tr>\n";

        
        print "<tr><td colspan=2>"
            . "<table width=100% >\n"
            . "<tr><th>variation</th><th>coordinates</th></tr>\n";
        foreach my $variant (@$variant_list_aref) {
            my ($type, $lend, $rend, $orient) = ($variant->{type},
                                                 $variant->{lend},
                                                 $variant->{rend},
                                                 $variant->{orient} );
            
            print "<tr><td>$type</td><td>$lend-$rend($orient)</td></tr>\n";
        }

        print "</table></td></tr>\n";
        print "</table>\n";

        &generate_3frame_image($reference_cdna_acc, $other_acc);

        
    }
}


####
sub generate_3frame_image {
    my ($reference_cdna_acc, $other_acc) = @_;
    
    my $ref_gene_obj = &get_gene_obj($reference_cdna_acc);
    my $other_gene_obj = &get_gene_obj($other_acc);

    ## Frame-based image
    print "<p>Illustration with exons tiered by genome sequence <b>frame</b></p>";
    my $gene_cdna_imager = new Gene_cdna_image;
    $gene_cdna_imager->{Sequence_feature_drawer_obj}->{DRAW_PANEL_SCALER} = 0.8;
    $gene_cdna_imager->{expand_3tiers} = "frame";
    
    my $file_token = "$reference_cdna_acc" . "_" . $other_acc;
    $file_token =~ s/\W/_/g;
    my $image_filename = "$ENV{WEBSERVER_TMP}/$file_token.frame.png";
    open (my $fh, ">$image_filename");
    binmode $fh;
    print $fh $gene_cdna_imager->create_image($ref_gene_obj, undef, $other_gene_obj)->png();
    close $fh;

    print "<img src=\"show_png.cgi?image=$image_filename\">\n\n";

=cds_phase_img_removed

    ## Phase-based image:
    print "<p>Illustration with exons tiered by incident cds <b>phase</b></p>";
    my $gene_cdna_imager = new Gene_cdna_image;
    $gene_cdna_imager->{Sequence_feature_drawer_obj}->{DRAW_PANEL_SCALER} = 0.8;
    $gene_cdna_imager->{expand_3tiers} = "phase";
    
    my $file_token = "$reference_cdna_acc" . "_" . $other_acc;
    $file_token =~ s/\W/_/g;
    my $image_filename = "$ENV{WEBSERVER_TMP}/$file_token.phase.png";
    open (my $fh, ">$image_filename");
    binmode $fh;
    print $fh $gene_cdna_imager->create_image($ref_gene_obj, undef, $other_gene_obj)->png();
    close $fh;

    print "<img src=\"show_png.cgi?image=$image_filename\">\n\n";

=cut
    
    print "<hr>\n";
    
}



####
sub get_gene_obj {
    my $acc = shift;

    my $gene_obj = &Ath1_cdnas::get_asmbl_gene_obj($dbproc, $acc, 0, 0);
    
    $gene_obj->{TU_feat_name} = $gene_obj->{Model_feat_name} = $acc;
        
    return ($gene_obj);
}
    
