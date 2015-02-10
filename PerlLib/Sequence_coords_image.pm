#!/usr/local/bin/perl

package Sequence_coords_image;

## Module developed by Brian Haas and Sam Angiuoli at The Institute for Genomic Research
## Simple API to generate images for displaying properties of biological sequences.
## bhaas@tigr.org, angiuoli@tigr.org

use GD;
use strict;
use vars qw ($DEBUG);

=NAME
    Sequence_coords_image.pm - A module to create images for the purpose of 
    illustrating the positions of any characteristic relative to a biological sequence.
    ie. hmm profiles, blast matches, promoter elements, binding sites, motifs, etc.

    A simple example program:
    
    # create a Sequence based image object
    my $obj = new Sequence_coords_image();
    # set the sequence length
    $obj->{SEQ_LENGTH} = 100;
    # create a single match element
    $element = new Sequence_coords_image::Element(25,"full", 50, "partial", "green");
    $obj->add_element($element);
    # create another match element
    $element = new Sequence_coords_image::Element(25,"full", 75, "partial", "blue");
    $obj->add_element($element);
    #generate and print the image.
    $img = $obj->create_image();
    print $img->png;

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
                     ELEMENTS => [], #array ref for holding each element (as an array reference, also);
                     SEQ_START => 0, #crop viewing area by setting start-stop
                     SEQ_STOP => 0, #if not specified use max coordinate
                     SEQ_LENGTH => 0, #DEPRECATED (left in for compatibility) use SEQ_START - SEQ_STOP instead
                     NO_BORDER => 0 #eliminates vertical white_space around image
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
sub add_element {
    my $self = shift;
    my $element_ref = shift;
    ## add params to the instance member ELEMENTS
    my $num_existing = $#{$self->{ELEMENTS}};
    $num_existing++;
    $self->{ELEMENTS}->[$num_existing] = $element_ref;
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
    my $num_elements = @{$self->{ELEMENTS}}; #+1
    if ($self->{TICKER_TOGGLE}) {
        $num_elements++;
    }
    
    if (!$self->{NO_BORDER}) {
        $num_elements++;
    }
    my $image_height = ($num_elements) * $element_vspacing; 
    ## create GD image.
    my $im = new GD::Image($image_x, $image_height);
    &set_initial_colors ($im,$self->{COLORS}); #populates color map #populates global %colors
    if($self->{SEQ_LENGTH}){
        $self->{SEQ_STOP} = $self->{SEQ_START}+$self->{SEQ_LENGTH};
    }
    else{
        unless ($self->{SEQ_STOP}) {
            $self->{SEQ_STOP} = &determine_max_length ($self);
        }
    }
    my $curr_y= $element_vspacing;
    if($self->{NO_BORDER}){
        $curr_y = 0;
    }
    if ($self->{TICKER_TOGGLE}) {
        $curr_y = &create_ticker ($self, $im);
    }
    if ($#{$self->{ELEMENTS}} >= 0) {
        &create_matchs_image($self, $im, $curr_y);
    }
    return ($im);
}


####
sub determine_max_length {
    my ($self) = @_;
    my @coordinates;
    my $elements_array_ref = $self->{'ELEMENTS'};
    foreach my $element_ref (@$elements_array_ref) {
        my @elements = @{$element_ref->{COORDS_TYPE_N_COLOR}};
        for (my $i = 0; $i <= $#elements; $i += 5) {
            push (@coordinates, $elements[$i], $elements[$i+2]);
        }
    }
    @coordinates = sort {$a<=>$b} @coordinates;
    return (pop(@coordinates));
}


####
sub create_matchs_image {
    my ($self, $im, $curr_y) = @_;
    my ($element_vspacing, $draw_panel_size, $seq_start,$seq_stop, $black) = ($self->{ELEMENT_VERTICAL_SPACING},
                                                                              $self->{DRAW_PANEL_SIZE},
                                                                              $self->{SEQ_START},
                                                                              $self->{SEQ_STOP},
                                                                              &get_colors("black", $im,$self->{COLORS}));
    my $TOTAL_VSPACE = int (0.5 * $element_vspacing);
    if($self->{NO_BORDER}==1){
        $TOTAL_VSPACE = $element_vspacing;
    }
    my $BEGIN_TEXT_X = $draw_panel_size; 
    foreach my $match_element (@{$self->{ELEMENTS}}) {
       	my $ycentral = int (0.5 * $TOTAL_VSPACE) + $curr_y;
        ## draw horizontal line thru center of the element position
        if ($match_element->{CENTRAL_LINE}) {
            $im->line (&coord_transform($self, [$seq_start, $ycentral, $seq_stop, $ycentral]), &get_colors ("black", $im,$self->{COLORS}));
        }
        my $tempref = $match_element->{COORDS_TYPE_N_COLOR};
        my @coordstuff = @$tempref;
        my($pos)=0;
        for (my $i = 0; $i <= $#coordstuff; $i+= 5) {
            my $end5 = $coordstuff[$i];
            my $end5_type = $coordstuff[$i+1];
            my $end3 = $coordstuff[$i+2];
            my $end3_type = $coordstuff[$i+3];
            my $color = $coordstuff[$i+4];
            
            #process arrowheads
            my(@img_coords) = &coord_transform($self, [$end5, $curr_y, $end3, ($curr_y+$TOTAL_VSPACE)]);
            my($cstr);
            foreach my $str (@img_coords){
                $cstr .= "$str,";
            }
            chop($cstr);
            $match_element->{IMAGE_COORDS}->{$end5} = $cstr;
            $pos++;
            if ($end5_type eq "arrow" || $end3_type eq "arrow") {
                &draw_arrow_element ($im, &coord_transform($self, [$end5, $curr_y, $end3, ($curr_y+$TOTAL_VSPACE)]), 
                                     &get_colors($color, $im,$self->{COLORS}), $element_vspacing, $end5_type, $end3_type);
            } else { #standard, just draw filled rectangle.
                $im->filledRectangle(&coord_transform($self,[$end5, $curr_y, $end3, 
                                                             ($curr_y + $TOTAL_VSPACE)]), &get_colors($color, $im,$self->{COLORS}));
            }
            
            # add broken glyph if partial
            if ($end5_type eq 'partial') {
                &draw_breakage ($im, $TOTAL_VSPACE, $draw_panel_size, &coord_transform($self, [$end5, $curr_y]),&get_colors("red", $im,$self->{COLORS}));
            }
            if ($end3_type eq 'partial') {
                &draw_breakage ($im, $TOTAL_VSPACE, $draw_panel_size, &coord_transform($self,[$end3, $curr_y]),&get_colors("red", $im,$self->{COLORS}));
            }
        }
        # add the within_match_text
        my @match_text = @{$match_element->{WITHIN_MATCH_TEXT}};
        while (@match_text) {
            my $coord = shift @match_text;
            my $text = shift @match_text;
            my $color = shift @match_text;
            $im->string(gdMediumBoldFont, &coord_transform($self, [$coord, $curr_y-$TOTAL_VSPACE]), $text, &get_colors($color, $im, $self->{COLORS}));
        }
        $im->string(gdSmallFont, $BEGIN_TEXT_X, $curr_y, $match_element->{TEXT}, $black);
        $curr_y += $element_vspacing;
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
    my ($self, $coords_ref) = @_;
    my (@coords) = @$coords_ref;
    my (@out_coords);
    my ($min_element_length,$max_element_length, $draw_panel_size, $draw_panel_spacing) = ($self->{SEQ_START},
                                                                                           $self->{SEQ_STOP} - $self->{SEQ_START}, 
                                                                                           $self->{DRAW_PANEL_SIZE}, 
                                                                                           $self->{DRAW_PANEL_SPACING});
    
    while (@coords) {
        my $xcoord = shift @coords;
        my $ycoord = shift @coords;
        print "xcoord_in\t$xcoord\tycoord_in\t$ycoord\n" if $DEBUG;
        $xcoord = (($xcoord-$min_element_length)/$max_element_length * ($draw_panel_size - $draw_panel_spacing)) ;
        $xcoord += $draw_panel_spacing / 2;
        $xcoord = int ($xcoord);
        print "xcoord_out\t$xcoord\tycoord_out\t$ycoord\n" if $DEBUG;
       	push (@out_coords, $xcoord, $ycoord);
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
    $colors->{"yellow"} = $im->colorAllocate(255,255,0);
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
    #print " ($min_element_length,$max_element_length, $element_vspacing) \n";
    
    my $curr_y = $element_vspacing;
    my $black = &get_colors('black', $im,$self->{COLORS});
    $im->line (&coord_transform($self, [$min_element_length, $curr_y, $max_element_length, $curr_y]), $black);
    my $ticker_height_small = int (0.20 * $element_vspacing);
    my $ticker_height_large = int (0.50 * $element_vspacing);
    my $line_text_pointer = int (0.3 * $element_vspacing);
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
            #if(length($i) >= 6){
            #	$label = $i/6;
            #	($flabel) = ($label =~ /(...)/);
            #	$flabel .= "M";
            #    }
            #    elsif(length($i) >= 4){
            #	$label = $i/4;
            #	($flabel) = ($label =~ /(...)/);
            #	$flabel .= "K";
            #    }
            #    else{
            #	$flabel = "$i";
            #    }
            $im->string(gdSmallFont, &coord_transform ($self, [$i, ($curr_y - $element_vspacing)]), "$flabel", $black);
            $im->line (&coord_transform ($self, [$i, $curr_y, $i, ($curr_y - $line_text_pointer)]), $black);
        }
        if ($line_length) {
            $im->line (&coord_transform ($self, [$i, $curr_y, $i, ($curr_y + $line_length)]), $black);
        }
    }
    return ($curr_y + $element_vspacing);
}




#################################################################################################################



package Sequence_coords_image::Element;

sub new {
    shift;
    my $self = { COORDS_TYPE_N_COLOR => [],
                 CENTRAL_LINE => 1,
                 TEXT => "",
                 WITHIN_MATCH_TEXT => [],
                 IMAGE_COORDS =>{}
             };
    #order of inp_params = ( (end5, end5_type, end3, end3_type, color) , ....)
    if (@_) {
        $self->{COORDS_TYPE_N_COLOR} = \@_;
    }
    bless ($self);
    return $self;
}

sub set_coords_type_n_color {
    my $self = shift;
    $self->{COORDS_TYPE_N_COLOR} = \@_;
}

sub set_central_line {
    my $self = shift;
    my $setting = shift;
    $self->{CENTRAL_LINE} = $setting;
}

sub set_text {
    #adds text to text panel
    my $self = shift;
    my $text = shift;
    $self->{TEXT} = $text;
}

sub get_coords {
    my($self,$end5) = @_;
    return $self->{IMAGE_COORDS}->{$end5};
}


sub add_match_text {
    #adds text within the drawing
    my $self = shift;
    #expected format ( (coord, text, color) , ... );
    push (@{$self->{WITHIN_MATCH_TEXT}}, @_);
}

1; #return true for use


