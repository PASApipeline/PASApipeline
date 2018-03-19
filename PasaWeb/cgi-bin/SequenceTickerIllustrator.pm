package SequenceTickerIllustrator;

use strict;
use warnings;
use Carp;
use GD::SVG;

####
sub new {
    my $packagename = shift;
    my $constructor_params = shift;

    # default settings overrided by constructor_params
    my $self = {
        major_unit => undef,
        minor_unit => undef,
        label => "major",
        draw_label_flag => 1,
        color => "black",
        line_top_offset => 0.2, #default central
        
        major_unit_height => 50,
        minor_unit_height => 20, # height in pixels
        
        
    };

    if (ref $constructor_params) {
        foreach my $param (keys %$constructor_params) {
            if (exists $self->{$param}) {
                $self->{$param} = $constructor_params->{$param};
            }
            else {
                confess "error, $param is not an attribute";
            }
        }
    }
    
    bless ($self, $packagename);
    return ($self);

}


####
sub draw {
    my $self = shift;
    my ($image, $start, $end, $canvas_position, $coord_converter, $colors_href) = @_;
    
    
    ## set step, major unit and minor unit unless already specified.
    
    ## major unit
    my $major_unit = $self->{major_unit};
    unless ($major_unit) {
        my $length = $start - $end + 1;
        my $order_10 = length($length) - 1;
        $order_10--;
        
        if ($order_10 >= 1) {
            $major_unit = 10 ** $order_10;
        }
        else {
            # default to 10
            $major_unit = 10;
        }
    }
    my $minor_unit = $self->{minor_unit};
    if ( (! $minor_unit) || $minor_unit > $major_unit) {
        
        $minor_unit = $major_unit / 10;
        if ($minor_unit < 1) {
            $minor_unit = 1;
        }
        
    }
    
    my $step = $minor_unit;

    unless ($start >= 0 && $end && $end > $start) {
        confess "Error, need start and end values and end must be greater than start: start: $start, end: $end";
    }

    ## unwrap the canvas rectangle:
    my ($x1, $y1, $x2, $y2) = ($canvas_position->{lend},
                               $canvas_position->{top},
                               $canvas_position->{rend},
                               $canvas_position->{bottom});
    
    unless ($x1 >= 0 && $y1 >= 0 && $x2 >= 0 && $y2 >= 0
            && $x2 >= $x1
            && $y2 >= $y1
            ) {
        confess "Error, rectangle is improper: ($x1, $y1) - ($x2, $y2)";
    }
    my $rect_height = $y2 - $y1;
    
    ## first, draw the central line:
    my $y_line = int ($rect_height * $self->{line_top_offset}) + $y1;
    
    my $ticker_color = $colors_href->{color};
    
    $image->line($x1, $y_line, $x2, $y_line, $ticker_color);
    
   
    my $major_unit_height = $self->{major_unit_height};
    my $minor_unit_height = $self->{minor_unit_height};
    
    my $tick_label_type = $self->{label};

    ## draw ticks and add labels:
    for (my $i = $start; $i <= $end; $i += $step) {
        if ($major_unit && ($i - $start) % $major_unit == 0) {
            ## draw minor tick:
            my $x_pos = $coord_converter->convert_coords($i);
            $image->line($x_pos, $y_line, $x_pos, $y_line + $major_unit_height, $ticker_color);
            
            if ($i && $tick_label_type =~ /major/i) {
                # draw label:
                $image->string(gdSmallFont, $x_pos, $y_line + $major_unit_height, "$i", $ticker_color) unless $i == $end;
            }
            
        } elsif ( $minor_unit && ($i - $start) % $minor_unit == 0) {
            ## draw minor tick:
            my $x_pos = $coord_converter->convert_coords($i);
            $image->line($x_pos, $y_line, $x_pos, $y_line + $minor_unit_height, $ticker_color);
            if ($i && $tick_label_type =~ /minor/i) {
                # draw label:
                $image->string(gdSmallFont, $x_pos, $y_line + $minor_unit_height, "$i", $ticker_color) unless $i == $end;
            }
        }
    }
    
}


1; #EOM
        
    
    
