#!/usr/local/bin/perl

package main;
our ($DEBUG);

package Sequence_feature_drawer;


use GD;
use strict;

=NAME
    Sequence_coords_image.pm - A module to create images for the purpose of 
    illustrating the positions of any characteristic relative to a biological sequence.
    ie. hmm profiles, blast matches, promoter elements, binding sites, motifs, etc.

    A simple example program:
    
    my $seq_feat_drawer = new Sequence_feature_drawer(
                                                      #default image length in pixels
                                                      IMAGE_X_SIZE => $IMAGE_X_SIZE || 1000,  
                                                      #percentage of image to draw matches, 
                                                      #rest for text-annotation of match
                                                      
                                                      DRAW_PANEL_SCALER => 0.6, 
                                                      #amount of vertical space consumed by each element
                                                      
                                                      ELEMENT_VERTICAL_SPACING => 12, 
                                                      
                                                      # turn on ticker
                                                      TICKER_TOGGLE => 1 );

    ## add a feature row to the image:
    my $curr_row = new Sequence_feature_drawer::Row(); 
    my $element = new Sequence_feature_drawer::Element(@{$feature->{element}});
    $curr_row->add_element($element);
    $obj->add_row($curr_row);
    my $img = $obj->create_image();





=cut 

## Package globals:
    #my %colors; #Can't do this as package global.  See bug fix above.
                ##found that if you keep using GD->color_allocate with the same color, 
                ##eventually, you lose the color info.
                ##so, now it is done only once per color and works fine.
    
## Constructor
    sub new {
        shift;
        my %input_params = @_;
        # set up defaults
        my $self = { IMAGE_X_SIZE => 600,  #default image length
                     DRAW_PANEL_SCALER => 0.60, #percentage of image to draw matches, rest for text-annotation of match
                     ELEMENT_VERTICAL_SPACING => 15, #amount of vertical space consumed by each match
                     TICKER_TOGGLE => 1, #toggle for displaying a ticker for protein length, default is on.
                     ROWS => [], #array ref for holding row elements.
                     SEQ_START => 0, #crop viewing area by setting start-stop
                     SEQ_STOP => 0, #if not specified use max coordinate
                     SEQ_LENGTH => 0, #DEPRECATED (left in for compatibility) use SEQ_START - SEQ_STOP instead
                     NO_BORDER => 0, #eliminates vertical white_space around image
                     ELEMENT_IMAGE_MAPPING => {},  # hashref->{feature_ID}->{seq_coord} = image_coord
                     TICKER_HEIGHT => 20,
                     
                     #gd image parameters:
                     background_color => "216:216:216",
                     
                     row_groupings => [],
                 }; 
# can manually set parameters
        foreach my $key (%input_params) {
            if (exists ($self->{$key})) {
                $self->{$key} = $input_params{$key};
            }
        }
        $self->{COLORS}={};
        bless($self);
        return $self;
    }

####
sub add_row {
    my $self = shift;
    my $element_ref = shift;
    ## add params to the instance member ELEMENTS
    my $num_existing = $#{$self->{ROWS}};
    $num_existing++;
    $self->{ROWS}->[$num_existing] = $element_ref;
}

sub add_row_grouping {
    my $self = shift;
    my @rows = @_;
    push (@{$self->{row_groupings}}, [@rows]);
}

sub get_row_groupings {
    my $self = shift;
    return (@{$self->{row_groupings}});
}

####
sub create_image {
    my $self = shift;
    my $image_x = $self->{IMAGE_X_SIZE};
    my $draw_panel_scaler = $self->{DRAW_PANEL_SCALER};
    my $element_vspacing = $self->{ELEMENT_VERTICAL_SPACING};
    my $draw_panel_size = $draw_panel_scaler * $image_x;
    $self->{DRAW_PANEL_SIZE} = $draw_panel_size;
    my $text_panel_size = $image_x - $draw_panel_size;
    $self->{TEXT_PANEL_SIZE} = $text_panel_size;
    my $draw_panel_spacing = 0.10 * $draw_panel_size;
    $self->{DRAW_PANEL_SPACING} = $draw_panel_spacing;
    my $num_elements = @{$self->{ROWS}};
    my $ticker_height = $self->{TICKER_HEIGHT};
    if (! $self->{TICKER_TOGGLE}) {
        $ticker_height = 0;
        $num_elements++;
    }
    if (!$self->{NO_BORDER}) {
        $num_elements++;
    }
    my $image_height = ($num_elements+1) * $element_vspacing + $ticker_height; 
    ## create GD image.
    my $im = new GD::Image($image_x, $image_height);
    &set_initial_colors ($im,$self->{COLORS}); #populates color map #populates global %colors
    unless ($self->{SEQ_STOP} && $self->{SEQ_START}) {
        my ($min_coord, $max_coord) = &determine_min_and_max_length ($self);
        unless ($self->{SEQ_START}) {
            $self->{SEQ_START} = $min_coord;
        }
        unless ($self->{SEQ_STOP}) {
            $self->{SEQ_STOP} = $max_coord;
        }
    }
    
    
    my $curr_y= $element_vspacing;
    if($self->{NO_BORDER}){
        $curr_y = 0;
    }
    if ($self->{TICKER_TOGGLE}) {
        $curr_y = &create_ticker ($self, $im);
    }
    if ($#{$self->{ROWS}} >= 0) {
        &create_matchs_image($self, $im, $curr_y);
    }
    return ($im);
}


####
sub determine_min_and_max_length {
    my ($self) = @_;
    my @coordinates;
    my $row_list_aref = $self->{ROWS};
    foreach my $row (@$row_list_aref) {
        unless (ref $row) { next;} #empty row
        my $elements_array_ref = $row;
        foreach my $element_ref ($row->get_elements()) {
            
            my @elements = @{$element_ref->{COORDS_TYPE_N_COLOR}};
            for (my $i = 0; $i <= $#elements; $i += 5) {
                push (@coordinates, $elements[$i], $elements[$i+2]);
            }
        }
    }
    @coordinates = sort {$a<=>$b} @coordinates;
    my $max_coord = pop @coordinates;
    my $min_coord = shift @coordinates;
    return ($min_coord, $max_coord);
}


####
sub create_matchs_image {
    my ($self, $im, $curr_y) = @_;

    # fill background color
    my $background = &get_colors($self->{background_color}, $im, $self->{COLORS});
    $im->fill(1,1,$background);
    
    my ($element_vspacing, $draw_panel_size, $seq_start,$seq_stop, $black) = ($self->{ELEMENT_VERTICAL_SPACING},
                                                                              $self->{DRAW_PANEL_SIZE},
                                                                              $self->{SEQ_START},
                                                                              $self->{SEQ_STOP},
                                                                              &get_colors("black", $im,$self->{COLORS}));
    print "Start: $seq_start, Stop: $seq_stop\n" if $DEBUG;
    
    my $TOTAL_VSPACE = int (0.5 * $element_vspacing);
    if($self->{NO_BORDER}==1){
        $TOTAL_VSPACE = $element_vspacing;
    }
    my $BEGIN_TEXT_X = $draw_panel_size; 
    my $rownum = 0;
    foreach my $row (@{$self->{ROWS}}) {
        $rownum++;
        
        print "Adding row ($rownum)\n" if $DEBUG;
        unless (ref $row) { 
            #empty row
            $curr_y += $element_vspacing;
            next; 
        } 
        
        foreach my $match_element ($row->get_elements()) {
            my $feature_ID = $match_element->{feature_ID};
            print "Adding match element.\n" if $DEBUG;
            my $ycentral = int (0.5 * $TOTAL_VSPACE) + $curr_y;
           
            ## Central Line processing:
            ## draw horizontal line thru center of the element position
            if (my $central_line = $match_element->{CENTRAL_LINE}) {
                if ($central_line == 1) { 
                    ## solid line connecting eles
                    my ($lend, $rend) = $match_element->get_element_span();
                    print "adding central line ($lend, $rend).\n" if $DEBUG;
                    $im->line (&coord_transform($self, [$lend, $ycentral, $rend, $ycentral]), 
                               &get_colors ("black", $im,$self->{COLORS}));
                }
                elsif ($central_line == 2) {
                    ## dashed line across entire line
                    $im->dashedLine (&coord_transform($self, [$seq_start, $ycentral, $seq_stop, $ycentral]), 
                               &get_colors ("120:120:00", $im,$self->{COLORS}));
                }
            }
            


            my $tempref = $match_element->{COORDS_TYPE_N_COLOR};
            my @coordstuff = @$tempref;
            
            for (my $i = 0; $i <= $#coordstuff; $i+= 5) {
                my $end5 = $coordstuff[$i];
                my $end5_type = $coordstuff[$i+1];
                my $end3 = $coordstuff[$i+2];
                my $end3_type = $coordstuff[$i+3];
                my $color = $coordstuff[$i+4];
                print "end5: $end5, end3: $end3\n" if $DEBUG;
                
                
                my(@img_coords) = &coord_transform($self, [$end5, $curr_y, $end3, ($curr_y+$TOTAL_VSPACE)]);
                if ($feature_ID) {
                    my ($y_low, $y_high) = sort {$a<=>$b} ($img_coords[1], $img_coords[3]);
                    $self->{ELEMENT_IMAGE_MAPPING}->{$feature_ID}->{y_low} = $y_low;
                    $self->{ELEMENT_IMAGE_MAPPING}->{$feature_ID}->{y_high} = $y_high;
                    $self->{ELEMENT_IMAGE_MAPPING}->{"~row:$rownum"}->{y_low} = $y_low;
                    $self->{ELEMENT_IMAGE_MAPPING}->{"~row:$rownum"}->{y_high} = $y_high;
                    
                    $self->{ELEMENT_IMAGE_MAPPING}->{$feature_ID}->{$end5} = $img_coords[0];
                    $self->{ELEMENT_IMAGE_MAPPING}->{$feature_ID}->{$end3} = $img_coords[2];
                }
                
                #process arrowheads
                if ($end5_type eq "arrow" || $end3_type eq "arrow") {
                    &draw_arrow_element ($im, [@img_coords], 
                                         &get_colors($color, $im,$self->{COLORS}), $element_vspacing, $end5_type, $end3_type);
                } else { #standard, just draw filled rectangle.
                    $im->filledRectangle(@img_coords, &get_colors($color, $im,$self->{COLORS}));
                }
                
                # add broken glyph if partial
                if ($end5_type eq 'partial') {
                    &draw_breakage ($im, $TOTAL_VSPACE, $draw_panel_size, &coord_transform($self, [$end5, $curr_y], $feature_ID),&get_colors("red", $im,$self->{COLORS}));
                }
                if ($end3_type eq 'partial') {
                    &draw_breakage ($im, $TOTAL_VSPACE, $draw_panel_size, &coord_transform($self,[$end3, $curr_y], $feature_ID),&get_colors("red", $im,$self->{COLORS}));
                }
            }
            # add the within_match_text
            my @match_text = @{$match_element->{WITHIN_MATCH_TEXT}};
            my $text_color = $match_element->get_text_color() || "black";
            
            while (@match_text) {
                my $coord = shift @match_text;
                my $text = shift @match_text;
                $im->string(gdSmallFont, &coord_transform($self, [$coord, $curr_y], $feature_ID), $text, &get_colors($text_color, $im, $self->{COLORS}));
            }
            $im->string(gdSmallFont, $BEGIN_TEXT_X, $curr_y, $match_element->{TEXT}, &get_colors($text_color, $im, $self->{COLORS}));
            
        }
        $curr_y += $element_vspacing;
    }
    if ($self->get_row_groupings()) {
        $self->connect_row_groupings($im);
    }
    

}


####
sub draw_arrow_element {
    my ($im, $x1, $y1, $x2, $y2, $color, $element_vspacing, $end5_type, $end3_type) = @_;
    my $max_arrowhead_height = int (0.8 * $element_vspacing);
    my $arrowhead_length = int (0.10 * ($x2 - $x1));
    if ($arrowhead_length < 2) {$arrowhead_length = int (0.7 *($x2 - $x1));} 
    my $arrowtip_rel_ycoord = int (($element_vspacing - $max_arrowhead_height) * .95) ;
    my $central_y = int ( ($y2 - $y1)/2 + $y1);
    my ($newx1, $newx2) = ($x1, $x2);
    if ($end3_type eq "arrow") {
        $newx2 = $x2 - $arrowhead_length;
        #draw arrowhead
        my $poly = new GD::Polygon;
        $poly->addPt($newx2, $y1 - $arrowtip_rel_ycoord);
        $poly->addPt($x2, $central_y);
        $poly->addPt($newx2, $y2 + $arrowtip_rel_ycoord);
        $im->filledPolygon($poly, $color);
    } 
    if ($end5_type eq "arrow") {
        $newx1 = $x1 + $arrowhead_length;
        #draw arrowhead
        my $poly = new GD::Polygon;
        $poly->addPt($newx1, $y1 - $arrowtip_rel_ycoord);
        $poly->addPt($x1, $central_y);
        $poly->addPt($newx1, $y2 + $arrowtip_rel_ycoord);
        $im->filledPolygon($poly, $color);
    }
    $im->filledRectangle($newx1, $y1, $newx2, $y2,$color);
}


####
sub draw_breakage {
    my ($im, $TOTAL_VSPACE, $draw_panel_size, $xcoord, $ycoord, $color) = @_;
    my $poly = new GD::Polygon;
    my $segmentV = int ($TOTAL_VSPACE/3);
    my $segmentH = int ($TOTAL_VSPACE * .2);
    $poly->addPt($xcoord, $ycoord);
    $poly->addPt($xcoord -= $segmentH, $ycoord += $segmentV);
    $poly->addPt($xcoord += $segmentH, $ycoord += $segmentV);
    $poly->addPt($xcoord -= $segmentH, $ycoord += $segmentV);
    $poly->addPt($xcoord += $segmentH, $ycoord);
    $poly->addPt($xcoord += $segmentH, $ycoord -= $segmentV);
    $poly->addPt($xcoord -= $segmentH, $ycoord -= $segmentV);
    $poly->addPt($xcoord += $segmentH, $ycoord -= $segmentV);
    $im->filledPolygon($poly, $color);#&get_colors("red", $im,$self->{COLORS}));
}



####
sub coord_transform {
    my ($self, $coords_ref, $feature_ID) = @_;
    my (@coords) = @$coords_ref;
    my (@out_coords);
    my ($min_element_length,$max_element_length, $draw_panel_size, $draw_panel_spacing) = ($self->{SEQ_START},
                                                                                           $self->{SEQ_STOP} - $self->{SEQ_START}, 
                                                                                           $self->{DRAW_PANEL_SIZE}, 
                                                                                           $self->{DRAW_PANEL_SPACING});
    
    while (@coords) {
        my $xcoord = shift @coords;
        my $ycoord = shift @coords;
        my $incoming_xcoord = $xcoord;
        #print "xcoord_in\t$xcoord\tycoord_in\t$ycoord\n" if $DEBUG;
        $xcoord = (($xcoord-$min_element_length)/$max_element_length * ($draw_panel_size - $draw_panel_spacing)) ;
        $xcoord += $draw_panel_spacing / 2;
        $xcoord = int ($xcoord);
        #print "xcoord_out\t$xcoord\tycoord_out\t$ycoord\n" if $DEBUG;
       	push (@out_coords, $xcoord, $ycoord);
        if ($feature_ID) {
            # store in mapping:
            $self->{ELEMENT_IMAGE_MAPPING}->{$feature_ID}->{$incoming_xcoord} = $xcoord;
        }
        
    }
    return (@out_coords);
} 


sub set_initial_colors {
    my ($im,$colors) = @_;
    $colors->{'white'} = $im->colorAllocate (255,255,255);
    #$im->transparent($colors->{'white'}); 
    $colors->{'black'} = $im->colorAllocate (0,0,0);
    $colors->{'red'} = $im->colorAllocate (255,0,0);
    $colors->{'blue'} = $im->colorAllocate (0,0,255);
    $colors->{'green'} = $im->colorAllocate (0,255,0);
}

####
sub get_colors {
    #manipulates global %colors
    my ($color, $im,$colors) = @_;
    my ($r, $g, $b);
    my $ret_color = $colors->{'black'};
    if (exists($colors->{$color})) {
        $ret_color = $colors->{$color};
    } else {
        if ($color =~ /^(\d+):(\d+):(\d+)/) {
            ($r, $g, $b) = ($1, $2, $3);
            if ( ($r>=0&&$r<=255) &&  ($g>=0&&$g<=255) &&  ($b>=0&&$b<=255)) {
                $colors->{$color} = $im->colorResolve($r, $g, $b);
                $ret_color = $colors->{$color};
            }
        }
    }
    return ($ret_color);
}

####
sub create_ticker {
    my ($self, $im) = @_;
    my ($min_element_length,$max_element_length, $element_vspacing) = ($self->{SEQ_START},$self->{SEQ_STOP}, $self->{ELEMENT_VERTICAL_SPACING});
    
    my $ticker_height = $self->{TICKER_HEIGHT};
    my $curr_y = $ticker_height;
    my $black = &get_colors('black', $im,$self->{COLORS});
    $im->line (&coord_transform($self, [$min_element_length, $curr_y, $max_element_length, $curr_y]), $black);
    my $ticker_height_small = int (0.20 * $ticker_height);
    my $ticker_height_large = int (0.50 * $ticker_height);
    my $line_text_pointer = int (0.3 * $ticker_height);
    my $value = int(($max_element_length-$min_element_length)/10);
    my $length = length($value);
    #print "value: $value\tLength: $length\n";
    my $lrg_interval = 10 ** ($length);
    #print "TICKER: $lrg_interval\n";
    if ($lrg_interval < 100) {$lrg_interval = 100;}
    my $sm_interval = int ($lrg_interval/10);
    my ($text_interval);
    if (int ((int ($max_element_length)-int($min_element_length))/$lrg_interval) <= 2) {
        $text_interval = 4 * $sm_interval;
    } else {
        $text_interval = $lrg_interval;
    }
    


    my($min_element_start) = $min_element_length - ($min_element_length % $sm_interval);
    my($max_element_start) = $max_element_length - ($max_element_length % $sm_interval) + $sm_interval;
    
    if (length($max_element_start) >= 6) {
        $text_interval = 2 * $lrg_interval;
    }
    

    for (my $i = $min_element_start; $i <= $max_element_start; $i+= $sm_interval) {
        #print "$i\n";
        my $line_length = 0;
        if ($i%$sm_interval == 0) {
            $line_length = $ticker_height_small;
        } 
        if ($i%$lrg_interval == 0) {
            $line_length = $ticker_height_large;
        }
        if ( $i%$text_interval == 0) {
            #add ticker text
            my($label,$flabel);
            $flabel = $i;
                        
            $im->string(gdSmallFont, &coord_transform ($self, [$i, ($curr_y - $ticker_height)]), "$flabel", $black);
            $im->line (&coord_transform ($self, [$i, $curr_y, $i, ($curr_y - $line_text_pointer)]), $black);
        }
        if ($line_length) {
            $im->line (&coord_transform ($self, [$i, $curr_y, $i, ($curr_y + $line_length)]), $black);
        }
    }
    return ($curr_y + 2*$element_vspacing);
}


####
sub connect_row_groupings {
    my $self = shift;
    my $im = shift;

    my @row_groupings = $self->get_row_groupings();
    
    foreach my $row_group (@row_groupings) {
        ## get coordinates and map to feat_ids
        my @eles;
        foreach my $row (@$row_group) {
            foreach my $match_element ($row->get_elements()) {
                my $feature_ID = $match_element->{feature_ID};
                my @coordstuff = @{$match_element->{COORDS_TYPE_N_COLOR}};
                for (my $i = 0; $i <= $#coordstuff; $i+= 5) {
                    my $end5 = $coordstuff[$i];
                    my $end3 = $coordstuff[$i+2];
        
                    ($end5, $end3) = sort {$a<=>$b} ($end5, $end3);
                    push (@eles, { end5 => $end5,
                                   end3 => $end3,
                                   feat_id => $feature_ID } );
                }
            }
        }

        ## link adjacent entries
        @eles = sort {$a->{end5}<=>$b->{end5}
                  ||
                  $a->{end3}<=>$b->{end3} 
                  } @eles;
        
        my $curr_ele = shift @eles;
        while (@eles) {
            my $next_ele = shift @eles;
            
            my $curr_feat_id = $curr_ele->{feat_id};
            my $curr_end5 = $curr_ele->{end5};
            my $curr_end3 = $curr_ele->{end3};

    
            my $next_feat_id = $next_ele->{feat_id};
            my $next_end5 = $next_ele->{end5};
            my $next_end3 = $next_ele->{end3};

            ## don't draw lines if coords overlap
            unless ($next_end5 < $curr_end3 && $next_end3 > $curr_end5) {
                #overlap
                
                
                ## get the coordinates:
                my $x1 = $self->{ELEMENT_IMAGE_MAPPING}->{$curr_feat_id}->{$curr_end3};
                my $y1_high = $self->{ELEMENT_IMAGE_MAPPING}->{$curr_feat_id}->{y_high};
                my $y1_low = $self->{ELEMENT_IMAGE_MAPPING}->{$curr_feat_id}->{y_low};
                my $y1 = int (($y1_high + $y1_low) /2);
                
                my $x2 = $self->{ELEMENT_IMAGE_MAPPING}->{$next_feat_id}->{$next_end5};
                my $y2_high = $self->{ELEMENT_IMAGE_MAPPING}->{$next_feat_id}->{y_high};
                my $y2_low = $self->{ELEMENT_IMAGE_MAPPING}->{$next_feat_id}->{y_low};
                my $y2 = int (($y2_high + $y2_low) /2);
                
                # draw connecting line
                $im->line($x1, $y1, $x2, $y2, &get_colors("black", $im, $self->{COLORS}));
            }
            
            # swap curr w/ next
            $curr_ele = $next_ele;
        }
        
    }
    

    
}

#################################################################################################################


package Sequence_feature_drawer::Row;
use strict;

sub new {
    my $packagename = shift;
    my $self = { ELEMENTS => [],
                 TEXT => "" };
    
    bless ($self, $packagename);
    return ($self);
}

sub add_element {
    my $self = shift;
    my @elements = @_;
    push (@{$self->{ELEMENTS}}, @elements);
}

sub get_elements {
    my $self = shift;
    return (@{$self->{ELEMENTS}});
}

sub set_text {
    my $self = shift;
    my $text = shift;
    $self->{TEXT} = $text;
}




######################################################################################################

package Sequence_feature_drawer::Element;

sub new {
    my $packagename = shift;
    my $self = { COORDS_TYPE_N_COLOR => [],
                 CENTRAL_LINE => 1,
                 WITHIN_MATCH_TEXT => [],
                 TEXT_COLOR => undef,
                 IMAGE_COORDS =>{},
                 LEND=>0,
                 REND=>0,
                 feature_ID=>undef   #used to get image-mapping coordinates for a feature. Must be assigned directly before image generation.
                 };
    bless ($self, $packagename);
    #order of inp_params = ( (end5, end5_type, end3, end3_type, color) , ....)
    if (@_) {
        $self->{COORDS_TYPE_N_COLOR} = \@_;
        $self->set_element_span();
    }
    
    return $self;
}


sub set_feature_ID {
    my $self = shift;
    my $feat_ID = shift;

    $self->{feature_ID} = $feat_ID;
}


sub set_coords_type_n_color {
    my $self = shift;
    $self->{COORDS_TYPE_N_COLOR} = \@_;
    $self->set_element_span();
}

sub set_central_line {
    my $self = shift;
    my $setting = shift;
    $self->{CENTRAL_LINE} = $setting;
}

sub set_text_color {
    my $self = shift;
    my $color = shift;
    $self->{TEXT_COLOR} = $color;
}

sub get_text_color {
    my $self = shift;
    my $color = $self->{TEXT_COLOR};
    return ($color);
}
    


sub set_text {
    #adds text to text panel
    my $self = shift;
    my $text = shift;
    $self->{TEXT} = $text;
}

sub set_element_span {
    my $self = shift;
    my @coords;
    my ($c1, $edge1, $c2, $edge2, $color, @coordstuff);
    @coordstuff =  @{$self->{COORDS_TYPE_N_COLOR}};
    while (@coordstuff) {
        ($c1, $edge1, $c2, $edge2, $color, @coordstuff) = @coordstuff;
        push (@coords, $c1, $c2);
    }
    @coords = sort {$a<=>$b} @coords;
    my $lend = shift @coords;
    my $rend = pop @coords;
    $self->{LEND} = $lend;
    $self->{REND} = $rend;
    return ($lend, $rend);
}

sub get_element_span {
    my $self = shift;
    my ($lend, $rend) = ($self->{LEND}, $self->{REND});
    unless ($lend && $rend) {
        return ($self->set_element_span());
    }
    return ($lend, $rend);
}


sub add_match_text {
    #adds text within the drawing
    my $self = shift;
    #expected format (coord, text);
    my $num_existing_elements = $#{$self->{WITHIN_MATCH_TEXT}};
    while (@_) {
        $num_existing_elements++;
        $self->{WITHIN_MATCH_TEXT}->[$num_existing_elements] = shift;
    }
}

1; #return true for use


