#!/usr/local/bin/perl

package Coordinate_draw_converter;

use strict;
use warnings;
use Carp;

sub new {
    my $packagename = shift;
    my ($image_bounds_aref, $feature_bounds_aref) = @_;
    
    my ($lend_image, $rend_image, 
        $min_left_coord, $max_right_coord) = (@$image_bounds_aref,
                                              @$feature_bounds_aref);
    ## make sure sorted
    ($lend_image, $rend_image) = sort {$a<=>$b} ($lend_image, $rend_image);
    ($min_left_coord, $max_right_coord) = sort {$a<=>$b} ($min_left_coord, $max_right_coord);

    my $self = {
        lend_image => $lend_image,
        rend_image => $rend_image,

        lend_coord => $min_left_coord,
        rend_coord => $max_right_coord,
    
    };
    
    bless ($self, $packagename);

    return ($self);
}
    

####
sub convert_coords {
    my $self = shift;
    my @coords = @_;
    
    my @retcoords;

    my ($lend_image, $rend_image) = ($self->{lend_image}, $self->{rend_image});
    my ($lend_coord, $rend_coord) = ($self->{lend_coord}, $self->{rend_coord});
    
    my $image_width = $rend_image - $lend_image + 1;
    my $coord_span = $rend_coord - $lend_coord + 1;

    foreach my $coord (@coords) {
        
        my $delta = $coord - $lend_coord;
        my $ratio = $delta / $coord_span;
        my $image_delta = $ratio * $image_width;
        my $image_coord = $lend_image + $image_delta;
        push (@retcoords, $image_coord);
    }
    
    
    if ( (! wantarray) && scalar (@retcoords) == 1) {
        return ($retcoords[0]);
    }
    else {
        return (@retcoords);
    }
}



1; #EOM
        
        
    
