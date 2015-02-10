package Splice_graph_illustrator;

use strict;
use warnings;
use Carp;
use Coordinate_draw_converter;

sub new {
    my $packagename = shift;

    my $params_href = shift;
    
    # some default settings that determine the positioning of the features.
    
    my $self = {
        feature_height => 20,
        arc_height => 30,
        arc_height_adjust => 20,
    };

    if (ref $params_href) {
        foreach my $param (keys %$params_href) {
            
            unless (exists $self->{$param}) {
                confess "Error, attribute $param is not an allowable attribute for obj.";
            }
            
            my $val = $params_href->{$param};
            $self->{$param} = $val;
        }
    }

    bless ($self, $packagename);

    return ($self);
}


####
sub draw_splice_graph {
    my $self = shift;
    my ($image, $nodes_aref, $coordinate_converter_obj, $bounding_box_href, $color_prefs_href) = @_;
    unless ($image 
            && ref $nodes_aref eq 'ARRAY' 
            && ref $coordinate_converter_obj 
            && ref $bounding_box_href eq 'HASH'
            && ref $color_prefs_href eq 'HASH') {
        confess "Error, params not properly set";
    }
    
    my $exon_color = $color_prefs_href->{exon} or confess "Need color for exon";
    my $intron_color = $color_prefs_href->{intron} or confess "Need color for intron";
    
    my ($box_lend, $box_rend, $box_top, $box_bottom) = ($bounding_box_href->{lend},
                                                        $bounding_box_href->{rend},
                                                        $bounding_box_href->{top},
                                                        $bounding_box_href->{bottom});
    
    ## separate into the introns and exons
    my @exon_nodes;
    my @intron_nodes;
        
    foreach my $node (@$nodes_aref) {
        my $type = $node->get_type();
        if ($type eq 'intron') {
            push (@intron_nodes, $node);
        }
        else {
            push (@exon_nodes, $node);
        }
    }
    
 
    my $midpt_Y = int ( ($box_top + $box_bottom) / 2);
    
    # get exon baselines
    my $exon_height = $self->{feature_height};
    my $exon_height_half = int($exon_height / 2);
    my $exon_top_Y = $midpt_Y - $exon_height_half;
    my $exon_bottom_Y = $midpt_Y + $exon_height_half;

    # draw the exons.
    foreach my $exon_node (@exon_nodes) {
        my ($exon_lend, $exon_rend) = $exon_node->get_coords();
        
        my ($exon_lend_coord, $exon_rend_coord) = $coordinate_converter_obj->convert_coords($exon_lend, $exon_rend);

        $image->filledRectangle($exon_lend_coord, $exon_top_Y, $exon_rend_coord, $exon_bottom_Y, $exon_color);

    }
    
    # draw the introns
    
    my $arc_height = $self->{arc_height};
    my $arc_height_adjust = $self->{arc_height_adjust};
    
    # alternate above and below
    my $pos = 1;
    foreach my $intron_node (sort {$a->{length}<=>$b->{length}}  @intron_nodes) {
        my ($intron_lend, $intron_rend) = $intron_node->get_coords();
        my ($intron_lend_coord, $intron_rend_coord) = $coordinate_converter_obj->convert_coords($intron_lend, $intron_rend);
        

        my $arc_midpt = int ( ($intron_lend_coord + $intron_rend_coord) / 2);
        my $arc_width = $intron_rend_coord - $intron_lend_coord;
        
        if ($pos) {
            # above
            $image->arc($arc_midpt, $exon_top_Y, $arc_width, $arc_height, 180, 0, $intron_color);
            
        }
        else {
            # below
            $image->arc($arc_midpt, $exon_bottom_Y, $arc_width, $arc_height, 0, 180, $intron_color);
            $arc_height += $arc_height_adjust;
        }
        $pos = ($pos == 1) ? 0 : 1; #alternate
        
    }
    
    

    return;
}
        
1; #EOM
