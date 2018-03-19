package main;

our ($IMAGE_X_SIZE, $DEBUG);

package Gene_cdna_image;
use strict;
use Sequence_feature_drawer;
use CDNA::CDNA_alignment;
use Gene_obj;


my $consensus = 20;
my $consensus_color = "0:172:187";

my $fl_color = "119:119:119";
my $non_fl_color = "32:33:48";

sub new {
    my $packagename = shift;

    my $self = { Sequence_feature_drawer_obj =>
                     new Sequence_feature_drawer(
                                                 IMAGE_X_SIZE => $IMAGE_X_SIZE || 750,  #default image length in pixels
                                                 DRAW_PANEL_SCALER => 0.6, #percentage of image to draw matches, rest for text-annotation of match
                                                 ELEMENT_VERTICAL_SPACING => 17, #amount of vertical space consumed by each element
                                                 TICKER_TOGGLE => 1 #toggle for displaying a ticker for protein length, default is on.
                                                 #SEQ_START => $min_end5
                                                 

                                                 ),
         
                                                 expand_3tiers => undef, # frame | phase    (frame is relative to genome (seqlen-1)%3,  phase is relative to complete CDS: the codon position of the first cds bp in the exon)  
                                                 
                                             };
    
    bless ($self, $packagename);

    return ($self);
}


sub create_image {
    my $self = shift;
    my (@alignmentsNGenes) = @_;
    my @data;
    foreach my $obj (@alignmentsNGenes) {
        my $struct;
        if ((ref $obj) =~ /Gene/) {
            $struct = $self->get_struct_from_gene_obj($obj);
        } elsif ((ref $obj) =~ /CDNA/) {
            $struct = $self->get_struct_from_alignment_obj($obj);
        }
        push (@data, $struct);
    }
    
    my %info;
    
    my @c;
    my $highlight = 0;
    foreach my $gene_struct (@data) {
        
        unless (ref $gene_struct) { next; } #may have empty entry for spacer

        my $coords_aref = $gene_struct->{coords};
        foreach my $href (@$coords_aref) {
            my $c1 = $href->{end5};
            my $c2 = $href->{end3};
            if (!$highlight) {
                $info{$c1} = 1;
                $info{$c2} = 1;
            }
            push (@c, $c1, $c2);
        }
        $highlight = 1; #only track coords for first entry.
    }
    @c = sort {$a<=>$b} @c;
    my $min_end5 = shift @c;
    $min_end5 -= 1;
    
    if ($DEBUG) {
        print "---\nMIN_COORD: $min_end5\n----\n";
        print "CREATING IMAGE\n";
    }
    
    
    my $obj = $self->{Sequence_feature_drawer_obj};  #store so can grab it from clients of module.
    
    my $SEPARATE_TIERS = $self->{expand_3tiers};
    
    foreach my $alignment (@data) {
        
        if (ref $alignment) { # otherwise, empty row.
            
            my @row_groupings;
            
            my $type = $alignment->{type};
            my @current_rows = (new Sequence_feature_drawer::Row());
            if ($SEPARATE_TIERS && $type eq "gene") {
                for (1..2) { ## add other 2 tiers
                    push (@current_rows, new Sequence_feature_drawer::Row());
                }
            }
                        
            my $accession = $alignment->{acc};
            my $name = $alignment->{name};
            
            my @coordsets = @{$alignment->{coords}};
                        
            
            my @lists = ([], [], []); #initialize 
            foreach my $coordset (@coordsets) {
                my @list;
                my $end5 = $coordset->{end5};
                my $end3 = $coordset->{end3};
                
                my $end5_agree = $info{$end5};
                my $end3_agree = $info{$end3};
                #$end5 -= $min_end5;
                #$end3 -= $min_end5;
                if ($DEBUG) {
                    print "END5: $end5\tEND3: $end3\t$accession\n";
                }
                my $color = ($accession =~ /\sFL/) ? $fl_color : $non_fl_color;
                
                push (@list, $end5, "full", $end3, "full", $color);
                if ((my $cds_end5 = $coordset->{cds_end5}) && (my $cds_end3 = $coordset->{cds_end3})) {
                    #$cds_end5 -= $min_end5;
                    #$cds_end3 -= $min_end5;
                    push (@list, $cds_end5, "full", $cds_end3, "full", "red");
                }
                if ($end5_agree) {
                    push (@list, $end5, "full", ($end5 + $consensus) , "full", $consensus_color);
                }
                if ($end3_agree) {
                    push (@list, ($end3 - $consensus), "full", $end3, "full", $consensus_color);
                }
                
                
                
                # element params (end5, end5_type, end3, end3_type, color, [... repeat for each match ....])
                #  end_types: (arrow, full, partial)
                # colors ( (white, blue, red, black, or green) or ("R:G:B") )
                if ($SEPARATE_TIERS && $type eq "gene") {
                    my $tier = $coordset->{tier};
                    unless (defined $tier) {
                        $tier = 1;
                    }
                    push (@{$lists[$tier]}, @list);
                }
                else {
                    push (@{$lists[0]}, @list);
                }
                
            }
            my $set_text_flag;
            for (my $i = 0; $i <= $#current_rows; $i++) {
                my $row = $current_rows[$i];
                my $element;
                my @list = @{$lists[$i]};
                #if (@list) {
                    $element = new Sequence_feature_drawer::Element(@list);
                #}
                if ($element) {
                    if ($SEPARATE_TIERS) {
                        $element->set_text($name . "-$SEPARATE_TIERS$i");
                        $element->set_feature_ID($accession . "-$SEPARATE_TIERS$i"); # add frame info
                        $element->{CENTRAL_LINE} = 2;
                    }
                    else {
                        $element->set_text($name);
                        $element->set_feature_ID($accession);
                    }
                    
                    if (my $polyAcoord = $alignment->{polyAcoord}) {
                        $element->add_match_text($polyAcoord, "A");
                        $element->set_text_color("green");
                    }
                    $row->add_element($element);
                }
                $obj->add_row($row);
                push (@row_groupings, $row);
            }
            $obj->add_row_grouping(@row_groupings);
        }
        else {
            # empty row
            $obj->add_row();
        }
    }
    
    #print image
    my $img = $obj->create_image();
    return ($img);
}


####
sub get_struct_from_gene_obj {
    my $self = shift;
    
    my $gene_obj = shift;
    my @coords;
    my $strand = $gene_obj->{strand};
    
    my @exons = $gene_obj->get_exons();
    my %cds;
    my $i = 0;
    
    foreach my $exon (@exons) {
        my ($end5, $end3) = sort {$a<=>$b} $exon->get_coords();
        $coords[$i]->{end5} = $end5;
        $coords[$i]->{end3} = $end3;
        if (my $cds_ref = $exon->get_CDS_obj()) {
            
            my ($end5, $end3) = sort {$a<=>$b} $cds_ref->get_coords();
            $coords[$i]->{cds_end5} = $end5;
            $coords[$i]->{cds_end3} = $end3;
            if (my $tier_type = $self->{expand_3tiers}) {
                if ($tier_type eq "frame") {
                    if ($strand eq '+') {
                        $coords[$i]->{tier} = ($end5 - 1 + 3 - $cds_ref->{phase}) % 3;
                    }
                    else { # (-) strand
                        $coords[$i]->{tier} = ($end5 - 1 + 3 + $cds_ref->{phase}) % 3;

                    }

                }
                elsif ($tier_type eq "phase") {
                    $coords[$i]->{tier} = $cds_ref->{phase};
                }
            }
        }
        $i++;
    }
    my $com_name = $gene_obj->{com_name} || "";
    my $struct = { type => 'gene',
                   strand => $strand,
                   coords => \@coords,
                   name => "($strand)" . $gene_obj->{Model_feat_name} . " $com_name",
                   acc => $gene_obj->{Model_feat_name},

                   };
    return ($struct);
    
}


####
sub get_struct_from_alignment_obj {
    my $self = shift;
   
    my $alignment = shift;
    my $acc = $alignment->get_acc();
    my $is_fli = $alignment->is_fli();
    my $fl_tag = ($is_fli) ? " FL" : "";
    my $strand = $alignment->get_orientation();
    my $spliced_orientation = $alignment->get_spliced_orientation();
    my @segments = $alignment->get_alignment_segments();
    my @coords;
    my $i=0;
    my $length = 0;
    foreach my $segment (@segments) {
        my ($end5, $end3) = $segment->get_coords();
        $coords[$i]->{end5} = $end5;
        $coords[$i]->{end3} = $end3;
        $i++;
        $length += abs($end3 - $end5) + 1;
    }
    my $type = "cdna";
    if ($acc =~ /^asmbl_/) {
        $type = "assembly";
        if ($fl_tag) {
            $fl_tag .= "-containing";
        }
    }
    my $title = $alignment->get_title() || "";
    
    my $struct =  { 
        type => $type,
        strand => $strand,
        coords => \@coords,
        name => "(a$strand/s$spliced_orientation) $acc $fl_tag $title",
        length=>$length,
        acc => $acc,
    };
    
    
    if (my $coord = $alignment->{polyAcoord}) {
        $struct->{polyAcoord} = $coord;
    }
    
    return ($struct);
}

1; #EOM
