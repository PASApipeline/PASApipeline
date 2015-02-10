#!/usr/bin/env perl

use GD;
use Data::Dumper;
use strict;
use Carp;
use lib ($ENV{EUK_MODULES});
use ColorGradient;

sub param {
    return (1);
}

our $DEBUG;

my $usage = "usage: $0 repeatFile [perID_filter start end]\n\n";
my $repeat_file = $ARGV[0] or die $usage;

my $start = $ARGV[2] || 0;
my $end   = $ARGV[3] || 0;
my $perid_filter = $ARGV[1] || 0;

if($start && $end) {
    ($start, $end) = sort {$a<=>$b} ($start, $end);
} 
else {
    ($start, $end) = (1,0);
}
my $error;

unless(-s $repeat_file) {
    $error = "Error: Can't find the repeat listings at: $repeat_file.<br>";
}

####################  
# Parse repeat file.
####################
open (REP, "$repeat_file") or croak "Error: Can't open $repeat_file\n";
my @repeats;
my $sequence_length;
my %seen;

while (<REP>) {
    unless (/\w/) {next;}
    chomp;
    my @x = split (/\s+/);
    my ($rep1_c1, $rep1_c2, $rep2_c1, $rep2_c2, $per_id, $seq_length) = 
		($x[1], $x[2], $x[5], $x[6], $x[8], $x[7]);
    
    ## check per_id, see if passes
    unless ($per_id >= $perid_filter) {next;} 
 
	## properly order:
	if ($rep2_c1 < $rep1_c1) {
		($rep1_c1, $rep1_c2, $rep2_c1, $rep2_c2) = ($rep2_c1, $rep2_c2, $rep1_c1, $rep1_c2); # swap 'em
	}
	
   
    #avoid dups.
    my @entries = sort (@x);
    my $key = join ("_", @entries);
    if ($seen{$key}) {
		next;
    } else {
		$seen{$key} = 1;
    }
    
    ## only examine repeats within specified limits.
    if ($end != 0) {
		unless (($rep1_c1 >= $start && $rep1_c1 <= $end) ||
				($rep2_c1 >= $start && $rep2_c1 <= $end) 
				) {
			next;
		}
    }
    
    my $rep = { rep1_c1 => $rep1_c1,
				rep1_c2 => $rep1_c2,
				rep2_c1 => $rep2_c1,
				rep2_c2 => $rep2_c2,
				per_id => $per_id};
    
	push (@repeats, $rep);
    $sequence_length = $seq_length;
}
close REP;

#####################
# Set sequence range.
#####################
unless ($end) {
    ### default setting.
    $end = $sequence_length;
}

#################
# Image drawing #
#################

## Set up the Scalable Image settings.
my $settings = { IMAGE_X_SIZE => 700,  #default image length
                 DRAW_PANEL_SCALER => 1.0, #percentage of image to draw matches, rest for text-annotation of match
                 ELEMENT_VERTICAL_SPACING => 10, #amount of vertical space consumed by each match
                 TICKER_TOGGLE => 1, #toggle for displaying a ticker for protein length, default is on.
				 ELEMENTS => [], #array ref for holding each element (as an array reference, also);
				 SEQ_START => $start, #crop viewing area by setting start-stop
				 SEQ_STOP => $end, #if not specified use max coordinate
				 DRAW_PANEL_SPACING => 50 #eliminates vertical white_space around image
				 }; 




my $image_x = $settings->{IMAGE_X_SIZE};
my $draw_panel_scaler = $settings->{DRAW_PANEL_SCALER};
my $element_vspacing = $settings->{ELEMENT_VERTICAL_SPACING};
my $draw_panel_size = $draw_panel_scaler * $image_x;
$settings->{DRAW_PANEL_SIZE} = $draw_panel_size;
my $text_panel_size = $image_x - $draw_panel_size;
my $image_height = 300;

## Calculate image_x coords for repeat features
foreach my $repeat (@repeats) {
    foreach my $type ("rep1_c1", "rep1_c2", "rep2_c1", "rep2_c2") {
		$repeat->{"${type}_tr"} = &coord_transform($settings, $repeat->{$type});
    }
}

## assign repeats to rows in the image:
my $max_row = &assign_repeats_to_rows(@repeats);
my ($im, $white, $blue, $black, $green, $red); #image and colors declared.
my $found_repeats = 0;
my $imagefile;
my $repeat_list_ref;
my $orient_list_ref;


 main: {
     ## create GD image.
     if (@repeats) {
		 $found_repeats = 1;
		 $im = new GD::Image($image_x, $image_height);
		 $white = $im->colorAllocate (255,255,255);
		 $blue = $im->colorAllocate (0,0,255);
		 $black = $im->colorAllocate (0,0,0);
		 $green = $im->colorAllocate (0,255,0);
		 $red = $im->colorAllocate (255,0,0);
		 
		 &create_ticker ($settings, $im);
		 my $section_num = 6;
		 &create_repeat_image($im, $section_num);
		 $section_num = 4; # 3 used up in ticker drawing.
		 
		 print $im->png();
     }
     
     exit(0);
 }     

sub create_repeat_image {
    my $im = shift;
    my $section_num = shift;
    my $baseline = shift;
    my $height = 0.75 * $settings->{ELEMENT_VERTICAL_SPACING};
    
    my $baseline = $image_height - ($section_num * $settings->{ELEMENT_VERTICAL_SPACING});
    
    #draw line where arcs will end.
    $im->line(&coord_transform($settings, $start), $baseline, &coord_transform($settings, $end), $baseline, $black);
    
    my $y_coord = sub { my $row_num = shift;
			return (int($image_height - ($section_num * $settings->{ELEMENT_VERTICAL_SPACING}) - ($row_num * $height)));
		    };
    
    my @arcs;
    my @rects;

	my @colors = &ColorGradient::get_RGB_gradient(scalar(@repeats));
	
    foreach my $repeat (@repeats) {
		my $color = shift @colors;
		my $draw_color = $im->colorAllocate(@$color);
		
		my $midpoint1 = ($repeat->{rep1_c1_tr} + $repeat->{rep2_c1_tr}) /2;
		my $width1 = abs (&coord_transform($settings, $repeat->{rep1_c1}) - &coord_transform($settings, $repeat->{rep2_c1})); 
		my $rownum1 = $repeat->{rep1_row};
		my $rep1_c1y = $y_coord->($rownum1);
		my $rep1_c2y = int($y_coord->($rownum1) + $height);
		$repeat->{rep1_c1y} = $rep1_c1y;
		$repeat->{rep1_c2y} = $rep1_c2y;

		my ($rep1_lend, $rep1_rend) = sort {$a<=>$b} ($repeat->{rep1_c1_tr},  $repeat->{rep1_c2_tr});
		
		print STDERR "BOX1: $rep1_c1y,$rep1_c2y\n";

		push (@rects, [($rep1_lend, $rep1_c1y, $rep1_rend, $rep1_c2y, $draw_color)]);
		push (@arcs, [$midpoint1,$baseline,
					  $width1,($baseline*2),
					  180,0,$draw_color]);
				
		my $midpoint2 = ($repeat->{rep1_c2_tr} + $repeat->{rep2_c2_tr}) /2;
		my $width2 = abs (&coord_transform($settings, $repeat->{rep1_c2}) - &coord_transform($settings, $repeat->{rep2_c2})); 
		my $rownum2 = $repeat->{rep2_row};
		my $rep2_c1y = $y_coord->($rownum2);
		my $rep2_c2y = int($y_coord->($rownum2) + $height);
		$repeat->{rep2_c1y} = $rep2_c1y;
		$repeat->{rep2_c2y} = $rep2_c2y;
		
		print STDERR "BOX2: $rep2_c1y,$rep2_c2y\n";
		
		my ($rep2_lend, $rep2_rend) = sort {$a<=>$b} ($repeat->{rep2_c1_tr}, $repeat->{rep2_c2_tr});
				
		push (@rects, [$rep2_lend, $rep2_c1y, $rep2_rend, 
					   $rep2_c2y, $draw_color]);
		
		push (@arcs, [$midpoint2,$baseline,
					  $width2,($baseline*2),
					  180,0,$draw_color]);
    }
    
    ## draw the arcs first,
    foreach my $arc (@arcs) {
		$im->arc(@$arc);
    }
    ## now draw the rectangles
    foreach my $rect (@rects) {
		$im->filledRectangle(@$rect);
		pop @$rect;
		push (@$rect, $black);
		$im->rectangle(@$rect);
    }
}



sub coord_transform {
    my ($settings, $xcoord, $pure_flag) = @_;
    my ($min_element_length,$max_element_length, $draw_panel_size, $draw_panel_spacing) = ($settings->{SEQ_START},
																						   $settings->{SEQ_STOP} - $settings->{SEQ_START}, 
																						   $settings->{DRAW_PANEL_SIZE}, 
																						   $settings->{DRAW_PANEL_SPACING});
    print "xcoord_in\t$xcoord\t" if $DEBUG;
    $xcoord = (($xcoord-$min_element_length)/$max_element_length * ($draw_panel_size - $draw_panel_spacing)) ;
    unless ($pure_flag) {
		$xcoord += $draw_panel_spacing / 2;
    }
    $xcoord = int ($xcoord);
    print "xcoord_out\t$xcoord\n" if $DEBUG;
    return ($xcoord);
} 


sub create_ticker {
    my ($settings, $im) = @_;
    my ($min_element_length,$max_element_length, $element_vspacing) = ($settings->{SEQ_START},$settings->{SEQ_STOP}, $settings->{ELEMENT_VERTICAL_SPACING});
    #print " ($min_element_length,$max_element_length, $element_vspacing) \n";
    my $curr_y = $image_height - 3 * $element_vspacing;
    $im->line (&coord_transform($settings, $min_element_length), $curr_y, &coord_transform($settings, $max_element_length), $curr_y, $black);
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
	    $im->string(gdSmallFont, &coord_transform ($settings, $i), ($curr_y + $element_vspacing), "$flabel", $black);
	    $im->line (&coord_transform ($settings, $i), $curr_y, &coord_transform($settings, $i), ($curr_y - $line_text_pointer), $black);
	}
	if ($line_length) {
	    $im->line (&coord_transform ($settings, $i), $curr_y, &coord_transform($settings, $i), ($curr_y + $line_length), $black);
	}
    }
    return ($curr_y + $element_vspacing);
}




sub assign_repeats_to_rows {
    my (@repeats) = @_;
    my @rows;
    my $max_row = 1;
    foreach my $repeat (@repeats) {
		foreach my $rep_type ("rep1", "rep2") {
			my (@coords) = sort {$a<=>$b} ($repeat->{"${rep_type}_c1"}, $repeat->{"${rep_type}_c2"});
			my $row_assignment;
			my $current_row = 1;
			while (!$row_assignment) {
				unless (ref $rows[$current_row]) {
					$rows[$current_row] = [];
				}
				my @currently_placed_elements = @{$rows[$current_row]};
				my $row_ok=1;
				foreach my $element (@currently_placed_elements) {
					my ($elem_c1, $elem_c2) = @$element;
					if ($coords[0] <= $elem_c2 && $coords[1] >= $elem_c1) { #overlap
						$row_ok = 0;
						last;
					}
				}
				if ($row_ok) {
					$row_assignment = $current_row;
					$repeat->{"${rep_type}_row"} = $row_assignment;
					push (@{$rows[$row_assignment]}, [@coords]);
				} else {
					$current_row++;
				}
			}
			if ($row_assignment > $max_row) { $max_row = $row_assignment;}
		}
    }
    return ($max_row);
}


