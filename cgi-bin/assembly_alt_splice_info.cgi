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
use CGI::Pretty ":standard";
use Pasa_CGI;


my $cgi = new CGI();
print $cgi->header('text/html');

my $db = $cgi->param('db');
my $cdna_acc = $cgi->param('cdna_acc');
my $DEBUG = $cgi->param('DEBUG');

unless ($cdna_acc && $db) {
    die "Need db, cdna_acc\n";
}


my $pid = $$;
unless ($ENV{WEBSERVER_TMP}) {
    $ENV{WEBSERVER_TMP} = "/tmp";
}

## image drawing settings:
my $image_x_size = 900;
my $draw_panel_scaler = 0.8;


my $show_all_flag = $cgi->param('SHOW_ALL');
my $show_alignments_flag = $cgi->param('SHOW_ALIGNMENTS');

my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $mysql_ro_user = &Pasa_conf::getParam("MYSQL_RO_USER");
my $mysql_ro_password = &Pasa_conf::getParam("MYSQL_RO_PASSWORD");

my ($dbproc) = &connect_to_db($mysql_server,$db,$mysql_ro_user,$mysql_ro_password);

my $css_common_text = &Pasa_CGI::get_common_CSS();


print $cgi->start_html(-title => "Alt-splicing report for $cdna_acc",
                       -head => style( { type => "text/css" }, $css_common_text ),
                       );


my $html = "";

my $query = "select annotdb_asmbl_id from clusters c, align_link al where c.cluster_id = al.cluster_id and al.align_acc = ?";

my $genomic_contig = &very_first_result_sql($dbproc, $query, $cdna_acc);

$html .= "<h1>Alternative Splicing Report for<br> [assembly: $cdna_acc, genome contig: $genomic_contig]</h1>\n";

$html .= "<h2>Splicing graph:</h2>\n";

## splicing graph image:
my $splicing_graph_img = "$ENV{WEBSERVER_TMP}/splice_graph.$pid.png";
$html .= "<img src=\"show_png.cgi?image=$splicing_graph_img\" alt=\'splicing graph for $cdna_acc\' >\n";


$html .= "<h2>Individual variations captured</h2>\n";
if ($show_all_flag) {
    $html .= "<font size=-2><a href=\"assembly_alt_splice_info.cgi?db=$db&cdna_acc=$cdna_acc&SHOW_ALIGNMENTS=$show_alignments_flag\">[hide underlying transcripts responsible]</a></font>\n";
} else {
    $html .= "<font size=-2><a href=\"assembly_alt_splice_info.cgi?db=$db&cdna_acc=$cdna_acc&SHOW_ALL=1&SHOW_ALIGNMENTS=$show_alignments_flag\">[show underlying transcripts responsible]</a></font>\n";
}
if ($show_alignments_flag) {
    $html .= "<font size=-2><a href=\"assembly_alt_splice_info.cgi?db=$db&cdna_acc=$cdna_acc&SHOW_ALL=$show_all_flag&SHOW_ALIGNMENTS=0\">[hide alignment text]</a></font>\n";
} else {
    $html .= "<font size=-2><a href=\"assembly_alt_splice_info.cgi?db=$db&cdna_acc=$cdna_acc&SHOW_ALL=$show_all_flag&SHOW_ALIGNMENTS=1\">[show alignment text]</a></font>\n";
}


## Retrieve variations and support
my $max_rend = 0;
my $got_alternate_internal_exons_flag = 0;
my @variations = &get_variations ($cdna_acc);  ## sets max_rend and got_alternate_exons_flag globals


if ($got_alternate_internal_exons_flag) {
    $html .= "<br><a href='#alternate_internal_exons'><font size=-1>view alternate internal exons report</font></a><br>\n";
}

my %cdna_acc_to_alignment_obj;
my %other_assembly_accs;

## generate variation summary table:
$html .= "<table border=1>\n";
$html .= "<tr><th>accession</th><th>splice variation</th><th>coordinates</th><th>inferred from</th></tr>\n";
my $count = 0;
foreach my $variation (sort {$a->{type} cmp $b->{type}} @variations) {
    my ($cdna_acc, $lend, $rend, $orient, $type, $support_list_aref) = ($variation->{cdna_acc},
                                                                        $variation->{lend},
                                                                        $variation->{rend},
                                                                        $variation->{orient},
                                                                        $variation->{type},
                                                                        $variation->{support_list}
                                                                        );
    
    $count++;
    my $image_filename = &create_alt_splice_image($variation, $count, $cdna_acc, @$support_list_aref);
    
    my $coordstring = "$lend-$rend($orient)";
    
            


    $html .= "<tr>"
        . "<td>$cdna_acc</td>"
        . "<td>$type</td>"
        . "<td>$coordstring</td>"
        . "<td>";
    foreach my $support_acc (@$support_list_aref) {
        $html .= "<a href=\"assembly_alt_splice_info.cgi?db=$db&cdna_acc=$support_acc\">$support_acc</a> ";
    }
    $html . "</td></tr>\n";
    
    $html .= "<tr><td colspan=4 ><img src=\"show_png.cgi?image=$image_filename\"></td></tr>\n";
    if ($show_alignments_flag) {
        ## add alignment text
        $html .= "<td colspan=4><font size=-2><table border=1>\n";
        foreach my $acc ($cdna_acc, @$support_list_aref) {
            my $alignment_text = &get_alignment_text($acc);
            
            $html .= "<tr><td>$acc $alignment_text</td></tr>\n";
        }
        $html .=  "</table></font></td>\n";
    }
        
    $html .= "<tr><td colspan=4>&nbsp;</td></tr>\n";

}

$html .=  "</table>\n";


&generate_splice_graph_image($splicing_graph_img, $cdna_acc, (keys %other_assembly_accs));




if ($got_alternate_internal_exons_flag) {
    $html .= &generate_alternate_internal_exons_report($cdna_acc);
}


print $html;
print $cgi->end_html();


exit(0);

####
sub create_alt_splice_image {
    my ($variation, $count, $cdna_acc, @other_accs) = @_;
    
    my $gene_cdna_imager = new Gene_cdna_image;

    $gene_cdna_imager->{Sequence_feature_drawer_obj}->{DRAW_PANEL_SCALER} = 0.8;
    $gene_cdna_imager->{Sequence_feature_drawer_obj}->{SEQ_STOP} = $max_rend;

    my @alignment_objs;

    foreach my $acc ($cdna_acc, undef, @other_accs) {
        
        
        if ($acc) {
            
            #print "<p>Getting alignment obj for $acc\n";
            
            
            my $alignment_obj = $cdna_acc_to_alignment_obj{$acc};
            unless (ref $alignment_obj) {
                
                #my $align_id = &Ath1_cdnas::get_validating_align_id_via_acc($dbproc, $acc);
                #$alignment_obj = $cdna_acc_to_alignment_obj{$acc} = 
                #    &Ath1_cdnas::create_alignment_obj($dbproc, $align_id);
                
                $alignment_obj = $cdna_acc_to_alignment_obj{$acc} = &Ath1_cdnas::get_alignment_obj_via_align_acc($dbproc, $acc);
                unless ($alignment_obj) {
                    die "Error, no alignment_obj for $acc";
                }
                
            }
            
            
            if ($acc =~ /^asmbl_/) { 
                ## pasa assembly
                my $gene_obj = &Ath1_cdnas::get_asmbl_gene_obj($dbproc, $acc, 1, 1);
                $gene_obj->{TU_feat_name} = $gene_obj->{Model_feat_name} = $acc;
                #print "<pre>" . $gene_obj->toString() . "</pre>\n";
                
                push (@alignment_objs, $gene_obj);
            }
            else {
                # transcript acc
                push (@alignment_objs, $alignment_obj);
            } 
            $other_assembly_accs{$acc} = 1;
        }
        else {
            #empty row
            push (@alignment_objs, undef);
        }
    }
    
    
    if ($show_all_flag && $variation) {
        my $sv_id = $variation->{sv_id};
        foreach my $other_acc (@other_accs) {
            my $query = "select transcripts_A, transcripts_B from splice_variation_support where sv_id = $sv_id and cdna_acc = \"$other_acc\"";
            #print "QUERY: $query\n";
            my $result = &first_result_sql($dbproc, $query);
            my ($transcripts_A, $transcripts_B) = @$result;
            
            push (@alignment_objs, undef); #empty spacer
            foreach my $acc (split (/,/, $transcripts_A), undef, split (/,/, $transcripts_B)) {
                if ($acc) {
                    #my $align_id = &Ath1_cdnas::get_validating_align_id_via_acc($dbproc, $acc);
                    #my $alignment_obj = $cdna_acc_to_alignment_obj{$acc} = 
                    #    &Ath1_cdnas::create_alignment_obj($dbproc, $align_id);
                    
                    my $alignment_obj = $cdna_acc_to_alignment_obj{$acc} = &Ath1_cdnas::get_alignment_obj_via_align_acc($dbproc, $acc);

                    push (@alignment_objs, $alignment_obj);
                } 
                else {
                    push (@alignment_objs, undef); #spacer
                }
            }
        }
    }
    
    
    my $seq_feature_drawer = $gene_cdna_imager->{Sequence_feature_drawer_obj};
    $seq_feature_drawer->{IMAGE_X_SIZE} = $image_x_size;
    $seq_feature_drawer->{DRAW_PANEL_SCALER} = $draw_panel_scaler;
    
    my $image = $gene_cdna_imager->create_image(@alignment_objs);

    my $manip_image = GD::Image->newFromPngData($image->png());
    
    if ($variation) {
        my ($lend, $rend, $orient, $type, $sv_id) = ( $variation->{lend},
                                                      $variation->{rend},
                                                      $variation->{orient},
                                                      $variation->{type},
                                                      $variation->{sv_id}
                                                      );
        
        
        
        
        my $y_high = $seq_feature_drawer->{ELEMENT_IMAGE_MAPPING}->{$cdna_acc}->{y_high};
        my $y_low = $seq_feature_drawer->{ELEMENT_IMAGE_MAPPING}->{$cdna_acc}->{y_low};
        
        my ($x1) = $seq_feature_drawer->coord_transform([$lend], $cdna_acc);
        my ($x2) = $seq_feature_drawer->coord_transform([$rend], $cdna_acc);
        
        my $purple = $manip_image->colorAllocate(129,0,183);
        $manip_image->rectangle($x1-2, $y_low-2, $x2+2, $y_high+2, $purple);
    }
    
    my $image_filename = "$ENV{WEBSERVER_TMP}/$pid.$count.$cdna_acc.alt_splice.png";
    open (my $fh, ">$image_filename") or die $!;
    binmode($fh);
    print $fh $manip_image->png();
    close $fh;
    
    return ($image_filename);
}


####
sub get_alignment_text {
    my $acc = shift;
    my $query = "select alignment from cdna_link where cdna_acc = ?";
    my $alignment = &very_first_result_sql($dbproc, $query, $acc);
    return ($alignment);
}



####
sub generate_splice_graph_image {
    my ($image_filename, @alignment_accs) = @_;  # alignment_acc includes targeted 'cdna_acc'

	unless (@alignment_accs) {
		die "Error, no alignment accs ";
	}

    my $element_vertical_spacing = 17;
    my $tier_vertical_spacing = 15;
    
    my $seq_feat_drawer = new Sequence_feature_drawer(
                                                      IMAGE_X_SIZE => $image_x_size, 
                                                      DRAW_PANEL_SCALER => $draw_panel_scaler,
                                                      ELEMENT_VERTICAL_SPACING => $element_vertical_spacing, 
                                                      TICKER_TOGGLE => 1,
                                                      SEQ_STOP => $max_rend,
                                                      );
    
    
    my %introns;
    
    
    my @alignment_objs;
    foreach my $other_acc (@alignment_accs) {
        my $other_alignment_obj = $cdna_acc_to_alignment_obj{$other_acc};
        unless (ref $other_alignment_obj) {
            $other_alignment_obj = $cdna_acc_to_alignment_obj{$other_acc} = &Ath1_cdnas::get_alignment_obj_via_align_acc($dbproc, $other_acc); 
            
        }
        push (@alignment_objs, $other_alignment_obj);
    }
    
    my $row = new  Sequence_feature_drawer::Row();

    my @list;
    
    my %exon_tokens;
    foreach my $alignment_obj (@alignment_objs) {
        my @segments = $alignment_obj->get_alignment_segments();
        foreach my $segment (@segments) {
            my ($lend, $rend) = $segment->get_coords();
            
            push (@list, $lend, "full", $rend, "full", "black");
            
            my $exon_token = "$lend" . "_" . "$rend";
            $exon_tokens{$exon_token} = 1;
        }
        
        my $element = new Sequence_feature_drawer::Element(@list);
        $element->set_feature_ID("splice_graph");
        $element->{CENTRAL_LINE} = 0;
        $row->add_element($element);
    }
    

    ## pull out introns for splice graph drawing:
    my %intron_coords;
    my @intron_features;
    

    my $target_alignment_obj = shift @alignment_objs or confess "Error, no alignemnt obj";
    
    my @target_intron_coords = $target_alignment_obj->get_intron_coords();
    foreach my $intron_coordset (@target_intron_coords) {
        my ($intron_lend, $intron_rend) = @$intron_coordset;
        my $intron_key = "$intron_lend" . "_" . "$intron_rend";
        $intron_coords{$intron_key} = "target";
        push (@intron_features, new TierFeatures::Feature($intron_lend, $intron_rend, $intron_key));
    }
    
    ## add other introns:
    my @other_intron_features;
    foreach my $other_alignment_obj (@alignment_objs) {
        my @other_intron_coords = $other_alignment_obj->get_intron_coords();
        foreach my $intron_coordset (@other_intron_coords) {
            my ($intron_lend, $intron_rend) = @$intron_coordset;
            my $intron_key = "$intron_lend" . "_" . "$intron_rend";
            unless ($intron_coords{$intron_key}) {
                # don't add it if it already exists:
                $intron_coords{$intron_key} = "other";
                push (@other_intron_features, new TierFeatures::Feature($intron_lend, $intron_rend, $intron_key));
            }
        }
    }
    
    # sort others by intron length
    @other_intron_features = sort { ($a->{rend} - $a->{lend}) <=> ($b->{rend} - $b->{lend}) } @other_intron_features;
    
    ## tier the intron features:
    my $feature_tierer = new TierFeatures();
    my @tiered_feats = $feature_tierer->tier_features(@intron_features, 
                                                      @other_intron_features);
    
    my $fh;
    if ($DEBUG) {
        open ($fh, ">/tmp/debug");
        print $fh Dumper (\@tiered_feats);
    }
    
    
    ## separate into top and bottom tiers
    my @top_tiers;
    my @bottom_tiers;
    my $count = 0;
    foreach my $tier (@tiered_feats) {
        $count++;
        if ($count % 2 == 0) {
            push (@bottom_tiers, $tier);
        } else {
            push (@top_tiers, $tier);
        }
    }

    @top_tiers = reverse @top_tiers;
    
    ## add empty row to for top tiers to draw in later
    foreach my $top_tier (@top_tiers) {
        $seq_feat_drawer->add_row();
    }
    
    ## add exon row
    $seq_feat_drawer->add_row($row);
    
    foreach my $bottom_tier (@bottom_tiers) {
        # add empty row 
        $seq_feat_drawer->add_row();
    }
    
    my $img = $seq_feat_drawer->create_image();
    my $manip_image = GD::Image->newFromPngData($img->png());
    
    ## draw the top intron layer:
    my $red = $manip_image->colorAllocate(255,0,0);
    my $green = $manip_image->colorAllocate(0,119,119);
    my $blue = $manip_image->colorAllocate(0,172,187);
    my $white = $manip_image->colorAllocate(255,255,255);
    
    
    my $target_intron_color = $green;
    my $other_intron_color = $red;
    
    

    my $num_top_tiers = scalar (@top_tiers);
    my $top_splice_graph_y = $seq_feat_drawer->{ELEMENT_IMAGE_MAPPING}->{"splice_graph"}->{y_low};
    my $bottom_splice_graph_y = $seq_feat_drawer->{ELEMENT_IMAGE_MAPPING}->{"splice_graph"}->{y_high};
    

    ## Draw the top tiers of introns:
    my $y_pos = $top_splice_graph_y - $num_top_tiers * $tier_vertical_spacing;


	#$manip_image->setThickness(1);

    foreach my $top_tier (@top_tiers) {
        my @intron_feats = @$top_tier;
        foreach my $intron_feat (@intron_feats) {
            my ($lend, $rend) = ($intron_feat->{lend}, $intron_feat->{rend});

            print $fh "drawing $lend-$rend\n" if $DEBUG;
            
            my $intron_key = "$lend" . "_" . "$rend";
            my $color = ($intron_coords{$intron_key} eq "target") 
                ? $target_intron_color
                : $other_intron_color;
            
            my ($x1) = $seq_feat_drawer->coord_transform([$lend], $intron_key);
            my ($x2) = $seq_feat_drawer->coord_transform([$rend], $intron_key);
            my $x_mid = int ( ($x1+$x2) / 2 + 0.5);
            
            $manip_image->line($x1, $top_splice_graph_y, $x_mid, $y_pos, $color);
            $manip_image->line($x2, $top_splice_graph_y, $x_mid, $y_pos, $color);

            ## draw ticks in exons:
            $manip_image->line($x1, $top_splice_graph_y, $x1, $bottom_splice_graph_y, $blue);
            $manip_image->line($x2, $top_splice_graph_y, $x2, $bottom_splice_graph_y, $blue);
            
        }
        $y_pos += $tier_vertical_spacing;
    }
    

    ## Draw the bottom tiers of introns:
    $y_pos = $bottom_splice_graph_y + $tier_vertical_spacing;
    
    foreach my $bottom_tier (@bottom_tiers) {
        my @intron_feats = @$bottom_tier;
        foreach my $intron_feat (@intron_feats) {
            my ($lend, $rend) = ($intron_feat->{lend}, $intron_feat->{rend});

            print $fh "drawing $lend-$rend\n" if $DEBUG;
            
            my $intron_key = "$lend" . "_" . "$rend";
            my $color = ($intron_coords{$intron_key} eq "target") 
                ? $target_intron_color
                : $other_intron_color;
            
            my ($x1) = $seq_feat_drawer->coord_transform([$lend], $intron_key);
            my ($x2) = $seq_feat_drawer->coord_transform([$rend], $intron_key);
            my $x_mid = int ( ($x1+$x2) / 2 + 0.5);
            
            $manip_image->line($x1, $bottom_splice_graph_y, $x_mid, $y_pos, $color);
            $manip_image->line($x2, $bottom_splice_graph_y, $x_mid, $y_pos, $color);
            
            ## draw ticks in exons:
            $manip_image->line($x1, $top_splice_graph_y, $x1, $bottom_splice_graph_y, $blue);
            $manip_image->line($x2, $top_splice_graph_y, $x2, $bottom_splice_graph_y, $blue);
            
        }
        $y_pos += $tier_vertical_spacing;
    }


    ## Draw white boxes around the exons:
    foreach my $exon_token (keys %exon_tokens) {
        my ($lend, $rend) = split (/_/, $exon_token);
        
        my ($x1) = $seq_feat_drawer->coord_transform([$lend], "splice_graph");
        my ($x2) = $seq_feat_drawer->coord_transform([$rend], "splice_graph");
        
        #$manip_image->rectangle($x1, $top_splice_graph_y, $x2, $bottom_splice_graph_y, $blue);
        $manip_image->line($x1, $top_splice_graph_y, $x2, $top_splice_graph_y, $blue);
        $manip_image->line($x1, $bottom_splice_graph_y, $x2, $bottom_splice_graph_y, $blue);
        
    }
    
    
    ## write the image file:
    open ($fh, ">$image_filename") or die $!;
    binmode $fh;
    print $fh $manip_image->png();
    close $fh;
}


####
sub generate_alternate_internal_exons_report {
    my $cdna_acc = shift;

    my $html_text = "<hr><a name='alternate_internal_exons'>&nbsp;</a>\n"
        . "<h2>Alternate Internal Exons Found</h2>\n";

    my @accs = ($cdna_acc);
    my $query = "select sv_id from splice_variation where cdna_acc = ? and type = 'retained_exon' and subtype = 'alternate_internal_exons' ";
    my @results = &do_sql_2D($dbproc, $query, $cdna_acc);
    
    my %other_accs;
    foreach my $result (@results) {
        my $sv_id = $result->[0];
        my $query = "select sv.cdna_acc from splice_variation sv, alt_splice_link asl where asl.sv_id_A = $sv_id and asl.sv_id_B = sv.sv_id";
        my $other_acc = &very_first_result_sql($dbproc, $query);
        unless ($other_acc) {
            print "Error, query $query didn't return entry!!!!\n";
            next;
        }
        $other_accs{$other_acc} = 1;
    }

    push (@accs, keys %other_accs);

    my $alt_internal_exon_image_name = "$ENV{WEBSERVER_TMP}/alt_internal_exons.$$.png";
    &generate_splice_graph_image($alt_internal_exon_image_name, @accs);

    $html_text .= "<p>Splice graph illustrating alternate internal exons:<br>\n";
    $html_text .= "<img src=\"show_png.cgi?image=$alt_internal_exon_image_name\">\n";
    
    shift @accs; # rid the cdna_acc entry
    
    ## draw separate images for each pair
    my $image_filename = &create_alt_splice_image (undef, ++$count, $cdna_acc, @accs);
    
    $html_text .= "<p>Relevant PASA assemblies:<br><img src=\"show_png.cgi?image=$image_filename\">\n";
    

    return ($html_text);

}


####
sub get_variations {
    my $cdna_acc = shift;
    
    my $query = "select sv_id, cdna_acc, lend, rend, orient, type, subtype from splice_variation where cdna_acc = ?";
    my @results = &do_sql_2D($dbproc, $query, $cdna_acc);
    
    foreach my $result (@results) {
        my ($sv_id, $cdna_acc, $lend, $rend, $orient, $type, $subtype) = @$result;
        my $support_list_aref = [];
        my $variation = { sv_id => $sv_id,
                          cdna_acc => $cdna_acc,
                          lend => $lend,
                          rend => $rend,
                          orient => $orient,
                          type => $type, 
                          subtype => $subtype,
                          support_list => $support_list_aref,
                          
                      };
        
        
        if ($subtype eq 'alternate_internal_exons') {
            $got_alternate_internal_exons_flag = 1;
        }
        
        my ($span_lend, $span_rend) = &Ath1_cdnas::get_alignment_span($dbproc, $cdna_acc);
        if ($span_rend > $max_rend) {
            $max_rend = $span_rend; # for image drawing.
        }
        
        my $query = "select cdna_acc from splice_variation_support where sv_id = ?";
        my @results = &do_sql_2D($dbproc, $query, $sv_id);
        foreach my $result (@results) {
            my ($other_acc) = @$result;
            push (@$support_list_aref, $other_acc);
            my ($span_lend, $span_rend) = &Ath1_cdnas::get_alignment_span($dbproc, $other_acc);
            
            if ($span_rend > $max_rend) {
                $max_rend = $span_rend; # for image drawing.
            }
            
        }
        
        push (@variations, $variation);
        
    }
    
    return (@variations);
}


    
    

    
    
