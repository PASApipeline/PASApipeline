package Gene_illustrator;

use strict;
use warnings;
use Carp;

sub new {
    my $packagename = shift;
    
    my $params_href = shift;
    
    # defaults
    my $self = {
        feature_height => 1.0, # ratio of total height.
    };
    
    if (ref $params_href) {
        foreach my $param (keys %$params_href) {
            
            unless (exists $self->{$param}) {
                confess "Error, attribute $param is not an allowable attribute for ob
j.";
            }
            
            my $val = $params_href->{$param};
            $self->{$param} = $val;
        }
    }
    
    bless ($self, $packagename);
    
    return ($self);
}


####
sub draw_gene {
    my $self = shift;
    my ($image, $gene, $coordinate_converter_obj, $bounding_box_href, $color_prefs_href) = @_;
    
    unless ($image 
            && ref $gene
            && ref $coordinate_converter_obj
            && ref $bounding_box_href eq 'HASH'
            && ref $color_prefs_href eq 'HASH') {
        confess "error, problem with params";
    }
    
    
    my ($box_lend, $box_rend, $box_top, $box_bottom) = ($bounding_box_href->{lend},
                                                        $bounding_box_href->{rend},
                                                        $bounding_box_href->{top},
                                                        $bounding_box_href->{bottom});
    
    my $box_height = $box_bottom - $box_top + 1;

    my $exon_color = $color_prefs_href->{exon};
    my $cds_color = $color_prefs_href->{CDS};

    my $feature_height_ratio = $self->{feature_height};

    my $offset = 0;
    if ($feature_height_ratio < 1) {
        $offset = int( ($box_height - ($feature_height_ratio * $box_height) ) / 2); 
    }

    my $top_Y = $box_top + $offset;
    my $bottom_Y = $box_bottom - $offset;
    
    my $midpt_Y = int ( ($bottom_Y + $top_Y) / 2);

    my ($gene_lend, $gene_rend) = $coordinate_converter_obj->convert_coords($gene->get_coords());

    # draw the central line
    $image->line($gene_lend, $midpt_Y, $gene_rend, $midpt_Y, $exon_color);

    my @exons = $gene->get_exons();
    foreach my $exon (@exons) {
        my ($exon_lend, $exon_rend) = sort {$a<=>$b} $coordinate_converter_obj->convert_coords($exon->get_coords());

        # draw filled rect
        $image->filledRectangle($exon_lend, $top_Y, $exon_rend, $bottom_Y, $exon_color);
        
        if (my $cds = $exon->get_CDS_exon_obj()) {
            my ($cds_lend, $cds_rend) = sort {$a<=>$b} $coordinate_converter_obj->convert_coords($cds->get_coords());

            $image->filledRectangle($cds_lend, $top_Y, $cds_rend, $bottom_Y, $cds_color);

        }
    }


    return;
}
        

1; #EOM
