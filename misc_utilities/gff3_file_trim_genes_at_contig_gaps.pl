#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Fasta_reader;
use GFF3_utils;
use Carp;
use Nuc_translator;
use Longest_orf;
use Data::Dumper;

my $usage = "\n\nusage: $0 gff3_file genome_db [genetic_code=universal|Tetrahymena|Candida]\n\n";

my $gff3_file = $ARGV[0] or die $usage;
my $fasta_db = $ARGV[1] or die $usage;
my $genetic_code = $ARGV[2];

my $MIN_GAP_LEN = 5;

#open (STDERR, ">&STDOUT");

my $write_gap_spanning_files = 0; #1;
my $VERBOSE = 0;


if ($genetic_code) {
    &Nuc_translator::use_specified_genetic_code($genetic_code);
}

my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
print STDERR "-loading GFF3 file: $gff3_file\n";
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

print STDERR "-loading genome from $fasta_db\n";
my $fasta_reader = new Fasta_reader($fasta_db);
my %genome = $fasta_reader->retrieve_all_seqs_hash();

my $TOTAL_MODELS_ADJUSTED = 0;

foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
    # print STDERR "// processing $asmbl_id\n";

    my $genome_seq = $genome{$asmbl_id} or die "Error, cannot locate sequence for genome entry: $asmbl_id";

    # print STDERR "-locating gaps in the sequence.\n";
    my @all_gaps = &find_contig_gaps($genome_seq);
    
    @all_gaps = sort {$a->[0]<=>$b->[0]} @all_gaps;
    
    my %gap_boundaries;
    foreach my $gap (@all_gaps) {
        my ($gap_lend, $gap_rend) = @$gap;
        $gap_boundaries{$gap_lend} = 1;
        $gap_boundaries{$gap_rend} = 1;
    }
    


    #print Dumper(\@all_gaps);
    
    my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    my @all_gene_objs;
    foreach my $gene_id (@gene_ids) {
        my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
        push (@all_gene_objs, $gene_obj_ref);
    }
    
    # print STDERR "-partitining genes based on exon overlap\n";
    my ($gene_objs_no_gap_exon_overlap_aref, $gene_objs_with_gap_exon_overlap_aref) = &partition_genes_based_on_gap_overlap(\@all_gaps, \@all_gene_objs);

    # unless (@$gene_objs_with_gap_exon_overlap_aref) { next; }
    
    print STDERR "There are " . scalar(@$gene_objs_no_gap_exon_overlap_aref) . " genes WITHOUT exon overlap\n";
    print STDERR "There are " . scalar(@$gene_objs_with_gap_exon_overlap_aref) . " genes WITH exon overlap\n";
    
    
    if (@$gene_objs_with_gap_exon_overlap_aref) {
        
        ## process the gap-overlapping gene models
        
        my @gap_overlapping_gene_models;
        foreach my $gene_obj_ref (@$gene_objs_with_gap_exon_overlap_aref) {
            
            ## KISS - remove isoforms for now, can regroup them later as needed.
            my @isoforms = $gene_obj_ref->get_additional_isoforms();
            
            $gene_obj_ref->delete_isoforms();
            $gene_obj_ref->refine_gene_object();
            $gene_obj_ref->set_CDS_phases(\$genome_seq);
            
            push (@gap_overlapping_gene_models, $gene_obj_ref);
            
            foreach my $isoform (@isoforms) {
                $isoform->refine_gene_object();
                $isoform->set_CDS_phases(\$genome_seq);
                
                push (@gap_overlapping_gene_models, $isoform);
            }
            
        }

        foreach my $gene_obj (@gap_overlapping_gene_models) {
            # I'm wishing this was already included in the gene obj
            my ($model_lend, $model_rend) = sort {$a<=>$b} $gene_obj->get_model_span();
            $gene_obj->{model_lend} = $model_lend;
            $gene_obj->{model_rend} = $model_rend;
        }

        if ($write_gap_spanning_files) {
            open (my $ofh, ">>gap-spanning_genes.$$.before.gff3") or die $!;
            foreach my $gene_obj (@gap_overlapping_gene_models) {
                print $ofh $gene_obj->to_GFF3_format() . "\n";
            }
            close $ofh;
        }
                
        ## 
        my @gene_objs = @gap_overlapping_gene_models;
        
        @gene_objs = sort {$a->{model_lend} <=> $b->{model_rend}} @gene_objs;
        
        my @genes_to_trim;
        my @genes_to_report;
        
        my %model_to_base_phases;
        
        while (@all_gaps) {
            my $gap = shift @all_gaps;
            my ($gap_lend, $gap_rend) = @$gap;
            
            print STDERR "GAP: $gap_lend-$gap_rend\n";
            
          gene_iter:
            while (@gene_objs) {
                if ($gene_objs[0]->{model_rend} < $gap_lend) {
                    # before current gap
                    print STDERR "Gene is before current gap: " . $gene_objs[0]->toString();
                    
                    my $gene_obj = shift @gene_objs;
                    push (@genes_to_report, $gene_obj);
                }
                elsif ($gene_objs[0]->{model_lend} <= $gap_rend) {
                    # must overlap current gap
                    my $gene_obj = shift @gene_objs;
                    push (@genes_to_trim, $gene_obj);
                }
                else {
                    last gene_iter;
                }
            }
            
            my @genes_to_trim_local = @genes_to_trim;
            @genes_to_trim = ();
            
            foreach my $gene_to_trim (@genes_to_trim_local) {
                if ($gene_to_trim->{model_rend} < $gap_lend) {
                    print STDERR "Gene to trim: " . $gene_to_trim->toString() . " is before gap: $gap_lend\n";
                    push (@genes_to_report, $gene_to_trim);
                }
                else {
                    my $model_feat_name = $gene_to_trim->{Model_feat_name};
                    $model_feat_name =~ s/----.*$//; # get original model name
                    unless (exists $model_to_base_phases{$model_feat_name}) {
                        $gene_to_trim->set_CDS_phases(\$genome_seq);
                        &assign_base_phase_encodings($model_feat_name, $gene_to_trim, \%model_to_base_phases);
                    }
                    
                    #print STDERR "Trimming gene at gap (" . join("-", @$gap) . "):\n" . $gene_to_trim->toString() . "\n\n";
                    my ($gene_pre_gap, $gene_post_gap) = &trim_gene_at_gap($gene_to_trim, $gap);
                    
                    if ($VERBOSE) {
                        print STDERR "Gene before gap: " . $gene_pre_gap->to_GFF3_format() . "\n" if $gene_pre_gap;
                        print STDERR "Gene post gap: " . $gene_post_gap->to_GFF3_format() . "\n" if $gene_post_gap;
                    }
                    
                    ## maybe just keep the longer of the two?
                    push (@genes_to_report, $gene_pre_gap) if $gene_pre_gap;
                    push (@genes_to_trim, $gene_post_gap) if $gene_post_gap; ## note, these are re-added to the queue for potential further trimming.
                    
                    $TOTAL_MODELS_ADJUSTED++;
                    
                }
            }
        }
        push (@genes_to_report, @genes_to_trim, @gene_objs);

        #print "Genes to report: @genes_to_report\n";
        #foreach my $gene_to_report (@genes_to_report) {
        #    print $gene_to_report->toString() . "\n";
        #}

        
        
        ## report them
        
        ## regroup by model.
        my %model_id_to_obj_list;
        foreach my $gene (@genes_to_report) {
            my $model_id = $gene->{Model_feat_name};
            $model_id =~ s/----.*$//;
            push (@{$model_id_to_obj_list{$model_id}}, $gene);
        }
        
        
        my @final_gene_objs;
        ## rejoin the trimmed parts and keep coding region in-frame
        foreach my $model_id (keys %model_id_to_obj_list) {
            
            my @gene_objs = @{$model_id_to_obj_list{$model_id}};
            if (scalar @gene_objs > 1) {
                print STDERR "// Got to trim:\n";
                foreach my $gene_obj (@gene_objs) {
                    print STDERR $gene_obj->toString();
                }
                print STDERR "\n\n";
                
                my $base_phase_info_href = $model_to_base_phases{$model_id} or die "Error, no base phase info for $model_id";
                
                my $reconstructed_gene = &merge_fragments($model_id, $base_phase_info_href, \@gene_objs, \%gap_boundaries);
                
                #print "Reconstructed:\n" . $reconstructed_gene->to_GFF3_format() . "\n";
                #die;
                push (@final_gene_objs, $reconstructed_gene);
            }
            else {
                push (@final_gene_objs, @gene_objs);
            }
            
            
        }
        
        ## Regroup by isoform.
        my %gene_id_to_obj_list;
        foreach my $gene_obj (@final_gene_objs) {
            
            my $gene_id = $gene_obj->{TU_feat_name};
            push (@{$gene_id_to_obj_list{$gene_id}}, $gene_obj);
            
        }
        
        
        foreach my $gene_id (keys %gene_id_to_obj_list) {
            my @isoforms = @{$gene_id_to_obj_list{$gene_id}};
            my $gene_obj = shift @isoforms;
            if (@isoforms) {
                $gene_obj->add_isoform(@isoforms);
                
                $gene_obj->refine_gene_object();
            }
            
            if (my $note = $gene_obj->{note}) {
                print "\n# $note\n";
            }
            
            print $gene_obj->to_GFF3_format() . "\n";
        
            if ($write_gap_spanning_files) {
                open (my $ofh, ">>gap-spanning_genes.$$.after.gff3") or die $!;
                
                if (my $note = $gene_obj->{note}) {
                    print "\n# $note\n";
                }
                
                print $ofh $gene_obj->to_GFF3_format() . "\n";
                close $ofh;
            }
            
        }
        
    }

    if (@$gene_objs_no_gap_exon_overlap_aref) {
        
        foreach my $gene_obj (@$gene_objs_no_gap_exon_overlap_aref) {
            print $gene_obj->to_GFF3_format() . "\n";
        }
    }
}

print STDERR "\n\n# Total models adjusted: $TOTAL_MODELS_ADJUSTED\n\n";

exit(0);


####
sub trim_gene_at_gap {
    my ($gene_to_trim, $gap) = @_;

    
    #print "//\n" . $gene_to_trim->toString() . "\nOverlaps gap: " . join("-", @$gap) . "\n";


    my ($gap_lend, $gap_rend) = @$gap;

    
    my @exon_coords;
    my @CDS_coords;

    my $orientation = $gene_to_trim->get_orientation();

    foreach my $exon ($gene_to_trim->get_exons()) {
        my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
        push (@exon_coords, [$lend, $rend]);
        
        if (my $cds_obj = $exon->get_CDS_obj()) {

            my ($cds_lend, $cds_rend) = sort {$a<=>$b} $cds_obj->get_coords();
            push (@CDS_coords, [$cds_lend, $cds_rend]);
        }
    }
    

    my ($before_trim_exon_coords_aref, $after_trim_exon_coords_aref) = &partition_coords_by_gap(\@exon_coords, [$gap_lend, $gap_rend]);
    my ($before_trim_cds_coords_aref, $after_trim_cds_coords_aref) = &partition_coords_by_gap(\@CDS_coords, [$gap_lend, $gap_rend]);

    my $before_gene_obj = &make_trimmed_gene_obj($gene_to_trim, $before_trim_exon_coords_aref, $before_trim_cds_coords_aref, "leftGap:$gap_lend-$gap_rend");
    my $after_gene_obj = &make_trimmed_gene_obj($gene_to_trim, $after_trim_exon_coords_aref, $after_trim_cds_coords_aref, "rightGap:$gap_lend-$gap_rend");
    

    #print "Left of gap: " . $before_gene_obj->toString() if $before_gene_obj;
    #print "Right of gap: " . $after_gene_obj->toString() if $after_gene_obj;
    
    return($before_gene_obj, $after_gene_obj);
    


}

sub find_contig_gaps {
    my ($sequence) = @_;
    
    my @gaps;

    while ($sequence =~ /(N+)/gi) {
        my $start = $-[0];
        my $stop = $+[0];
        
        $start++;
        
        my $gap_len = $stop - $start + 1;
        if ($gap_len >= $MIN_GAP_LEN) {
            push (@gaps, [$start, $stop] );
        }
    }

    return(@gaps);
}


####
sub partition_coords_by_gap {
    my ($exon_coords_aref, $gap_coords_aref) = @_;
    
    my ($gap_lend, $gap_rend) = @$gap_coords_aref;
    
    my @before_exon_coords;
    my @after_exon_coords;
    
    foreach my $exon_coordset (@$exon_coords_aref) {
        my ($exon_lend, $exon_rend) = @$exon_coordset;
        
        if ($exon_lend <= $gap_rend && $exon_rend >= $gap_lend) {
            my ($before_coords, $after_coords) = &trim_coords_to_gap([$exon_lend, $exon_rend], [$gap_lend, $gap_rend]);
            
            if (@$before_coords) {
                push (@before_exon_coords, $before_coords);
            }
            if (@$after_coords) {
                push (@after_exon_coords, $after_coords);
            }
        }
        elsif ($exon_rend < $gap_lend) {
            push (@before_exon_coords, $exon_coordset);
        }
        elsif ($exon_lend > $gap_rend) {
            push (@after_exon_coords, $exon_coordset);
        }
        else {
            die "Wassup!";
        }
        
    }
    
    ## require at least one codon in the trimmed region near the gap
    my @left_coords = @before_exon_coords;
    my @right_coords = @after_exon_coords;
    
    @left_coords = reverse sort {$a->[0]<=>$b->[0]} @left_coords;
    while (@left_coords) {
        my $terminal_seg_aref = $left_coords[0];
        my $len = abs($terminal_seg_aref->[1] - $terminal_seg_aref->[0]) + 1;
        if ($len < 3) {
            shift @left_coords;
        }
        else {
            last;
        }
    }
    
    @right_coords = sort {$a->[0]<=>$b->[0]} @right_coords;
    while (@right_coords) {
     my $terminal_seg_aref = $right_coords[0];
        my $len = abs($terminal_seg_aref->[1] - $terminal_seg_aref->[0]) + 1;
        if ($len < 3) {
            shift @right_coords;
        }
        else {
            last;
        }
    }   
    
    
    return(\@left_coords, \@right_coords);
}

####
sub trim_coords_to_gap {
    my ($exon_coords_aref, $gap_coords_aref) = @_;  # single pair of coordinates for each

    my ($exon_lend, $exon_rend) = @$exon_coords_aref;
    my ($gap_lend, $gap_rend) = @$gap_coords_aref;
    
    my @left_coords;
    my @right_coords;

    if ($exon_lend < $gap_lend) {
        my $left_coord_A = $exon_lend;
        my $left_coord_B = ($exon_rend < $gap_lend) ? $exon_rend : $gap_lend - 1;
        
        
        @left_coords = ($left_coord_A, $left_coord_B);
    }
    if ($exon_rend > $gap_rend) {
        my $right_coord_A = ($exon_lend > $gap_rend) ? $exon_lend : $gap_rend + 1;
        my $right_coord_B = $exon_rend;
        @right_coords = ($right_coord_A, $right_coord_B);
    }
    
    return(\@left_coords, \@right_coords);
}


####
sub make_trimmed_gene_obj {
    my ($initial_gene, $exon_coords_aref, $cds_coords_aref, $name_token) = @_;

    my $orientation = $initial_gene->get_orientation();

    
    my $gene_copy = $initial_gene->clone_gene();

    $gene_copy->erase_gene_structure();
    $gene_copy->{TU_feat_name} .= "----$name_token";
    $gene_copy->{Model_feat_name} .= "----$name_token";
    
    my %exon_coords;
    my %cds_coords;
    
    foreach my $exon_coordset (@$exon_coords_aref) {
        my ($end5, $end3) = @$exon_coordset;
        if ($orientation eq '-') {
            ($end5, $end3) = ($end3, $end5);
        }
        $exon_coords{$end5} = $end3;
    }
    
    foreach my $cds_coordset (@$cds_coords_aref) {
        my ($end5, $end3) = @$cds_coordset;
        if ($orientation eq '-') {
            ($end5, $end3) = ($end3, $end5);
        }
        $cds_coords{$end5} = $end3;
    }

    if (%cds_coords) {
        $gene_copy->populate_gene_object(\%cds_coords, \%exon_coords);
        
        my ($model_lend, $model_rend) = sort {$a<=>$b} $gene_copy->get_model_span();
        $gene_copy->{model_lend} = $model_lend;
        $gene_copy->{model_rend} = $model_rend;
        
        return($gene_copy);
    }
    else {
        return (undef);
    }
}


####
sub set_base_level_phase_info {
    my ($isoform) = @_;

    

    my @exons = $isoform->get_exons();

    foreach my $exon (@exons) {
        if (my $cds_exon_obj = $exon->get_CDS_obj()) {
            
            my ($end5, $end3) = $cds_exon_obj->get_coords();
            
            my $phase = $cds_exon_obj->{phase};
            print STDERR join("\t", $isoform->{Model_feat_name}, "$end5-$end3", $phase);
        }
    }

    return;
}

####
sub assign_base_phase_encodings {
    my ($model_feat_name, $gene_obj, $model_to_base_phases_href) = @_;
    
    print STDERR "ASSIGING BASE PHASE ENCODINGS FOR: $model_feat_name\n";
    

    my $orient = $gene_obj->get_orientation();
    my ($left_span, $right_span) = sort {$a<=>$b} $gene_obj->get_model_span();
    
    my @base_phases; # starting at left_span coordinate

    ## assign base level phase info
    foreach my $exon ($gene_obj->get_exons()) {
        
        if (my $cds_obj = $exon->get_CDS_exon_obj()) {

            my $local_cds_len = 0;
            my ($cds_end5, $cds_end3) = $cds_obj->get_coords();
            #print "CDS: $cds_end5-$cds_end3\n";
            my $cds_start_phase = $cds_obj->{phase};
            unless (defined $cds_start_phase) {
                die "Error, no CDS start phase defined for cds $cds_end5-$cds_end3 of gene: " . $gene_obj->toString();
            }
            
            $local_cds_len += $cds_start_phase;
            
            if ($orient eq '+') {
                for (my $i = $cds_end5; $i <= $cds_end3; $i++) {
                    my $phase = $local_cds_len % 3;
                    $base_phases[$i-$left_span] = $phase;
                    $local_cds_len++;
                }
            }
            elsif ($orient eq '-') {
                for (my $i = $cds_end5; $i >= $cds_end3; $i--) {
                    my $phase = $local_cds_len % 3;
                    $base_phases[$i-$left_span] = $phase;
                    $local_cds_len++;
                }
            }
        }
    }

    
    #foreach my $b (@base_phases) {
    #    unless (defined $b) {
    #        $b = '-';
    #    }
    #}


    #print join(", ", @base_phases) . "\n";
   

    $model_to_base_phases_href->{$model_feat_name} = { start_pos => $left_span,
                                                       base_phases => \@base_phases,
    };
    
    return;
}


####
sub merge_fragments {
    my ($model_id, $base_phase_info_href, $gene_objs_aref, $gap_boundaries_href) = @_;

    my $start_pos = $base_phase_info_href->{start_pos};
    my $base_phases_aref = $base_phase_info_href->{base_phases};
    
    my @gene_objs = @$gene_objs_aref;

    @gene_objs = sort {$a->{model_span}->[0] <=> $b->{model_span}->[0]} @gene_objs;

    my $orient = $gene_objs[0]->get_orientation();

    my $left_gene = shift @gene_objs;

    foreach my $right_gene (@gene_objs) {
        
        my $left_terminal_exon = &get_rightmost_terminal_exon($left_gene);
        my $right_terminal_exon = &get_leftmost_terminal_exon($right_gene);
        
        my $left_cds_obj = $left_terminal_exon->get_CDS_exon_obj() or die "Error, no CDS obj for right terminal exon of " . $left_gene->toString();
        my $right_cds_obj = $right_terminal_exon->get_CDS_exon_obj() or die "Error, no CDS obj for left terminal exon of " . $right_gene->toString();

        if ($orient eq '+') {
            
            my $left_cds_end3 = $left_cds_obj->{end3};
            my $right_cds_end5 = $right_terminal_exon->{end5};

            my $phase_left = $base_phases_aref->[ $left_cds_end3 - $start_pos ];
            my $phase_right = $base_phases_aref->[ $right_cds_end5 - $start_pos ];

            print STDERR "Phases to join (+): $phase_left, $phase_right\n";
            
            if ( ($phase_left + 1) % 3 == $phase_right) {
                ## OK, simple merge, no adjustment necessary

                print STDERR "Simple merging. No adjustment necessary.\n";
                $left_gene = &merge_genes($left_gene, $right_gene);

            }
            else {
                ## must adjust boundaries
                if ( $gap_boundaries_href->{$left_cds_end3 + 1} ) {
                    ## adjust left segment
                    print STDERR "** adjusting left segment (+).\n";
                    
                    $left_cds_end3--;
                    $left_terminal_exon->{end3}--;
                    $left_cds_obj->{end3}--;
                    
                    my $counter = 0;
                    while ( ($base_phases_aref->[$left_cds_end3 - $start_pos] + 1) % 3 != $phase_right) {
                        $left_cds_end3--;
                        $left_terminal_exon->{end3}--;
                        $left_cds_obj->{end3}--;
                        
                        $counter++;
                        if ($counter > 2) { die "too many edits";}
                                                
                    }
                
                }
                elsif ($gap_boundaries_href->{$right_cds_end5 - 1}) {
                    ## adjust right segment
                    print STDERR "** adjusting right segment (+).\n";
                    $right_cds_end5++;
                    $right_terminal_exon->{end5}++;
                    $right_cds_obj->{end5}++;
                    
                    my $counter = 0;
                    while ( ($phase_left + 1) % 3 != $base_phases_aref->[$right_cds_end5 - $start_pos]) {
                        $right_cds_end5++;
                        $right_terminal_exon->{end5}++;
                        $right_cds_obj->{end5}++;

                        $counter++;
                        if ($counter > 2) { die "too many edits";}
                    }


                }
                else {
                    print STDERR "Error, cannot figure out which segment neighbors the gap (+).";
                    $left_gene->{note} = "ERROR WITH MERGE";
                }
                
                $left_gene = &merge_genes($left_gene, $right_gene);
            }

            
        }
        else {
            # orient '-'
                        
            my $left_cds_end5 = $left_cds_obj->{end5};
            my $right_cds_end3 = $right_terminal_exon->{end3};

            my $phase_left = $base_phases_aref->[ $left_cds_end5 - $start_pos ];
            my $phase_right = $base_phases_aref->[ $right_cds_end3 - $start_pos ];

            print STDERR "Phases to join (-): $phase_right, $phase_left\n";
            
            if ( ($phase_right + 1) % 3 == $phase_left) {
                ## OK, simple merge, no adjustment necessary
                
                $left_gene = &merge_genes($left_gene, $right_gene);
                
            }
            else {
                ## must adjust boundaries
                if ($gap_boundaries_href->{$left_cds_end5 + 1}) {
                    print STDERR "** adjusting left segment (-).\n";
                    
                    $left_cds_end5--;
                    $left_terminal_exon->{end5}--;
                    $left_cds_obj->{end5}--;
                    
                    my $counter = 0;
                    while ($base_phases_aref->[$left_cds_end5 - $start_pos]  != ($phase_right + 1) % 3) {
                        $left_cds_end5--;
                        $left_terminal_exon->{end5}--;
                        $left_cds_obj->{end5}--;
                        
                        $counter++;
                        if ($counter > 2) { die "too many edits";}
                    }
                }
                elsif ($gap_boundaries_href->{$right_cds_end3 -1}) {
                    print STDERR "** adjusting right segment (-).\n";
                    
                    
                    $right_cds_end3++;
                    $right_terminal_exon->{end3}++;
                    $right_cds_obj->{end3}++;
                    
                    my $counter = 0;
                    while ($phase_left != ($base_phases_aref->[$right_cds_end3 - $start_pos] + 1) % 3 ) { 
                        $right_cds_end3++;
                        $right_terminal_exon->{end3}++;
                        $right_cds_obj->{end3}++;
                        
                        $counter++;
                        if ($counter > 2) { die "too many edits";}
                        
                    }
                    
                }
                else {
                    print STDERR "Error, cannot figure out which segment neighbors the gap (-).";
                    $left_gene->{note} = "ERROR WITH MERGE";
                }

                $left_gene = &merge_genes($left_gene, $right_gene);
            
            }
            

        }
        
        
        
    }
    
    return($left_gene);
}


####
sub get_rightmost_terminal_exon {
    my ($gene) = @_;

    my @exons = $gene->get_exons();
    @exons = sort {$a->{end5}<=>$b->{end5}} @exons;

    my $rightmost_exon = pop @exons;

    return($rightmost_exon);
}

####
sub get_leftmost_terminal_exon {
    my ($gene) = @_;

    my @exons = $gene->get_exons();
    @exons = sort {$a->{end5}<=>$b->{end5}} @exons;

    my $leftmost_exon = shift @exons;
    
    return($leftmost_exon);

}


####
sub merge_genes {
    my ($left_gene, $right_gene) = @_;

    ## just add the right gene's exons to the left gene

    my @right_gene_exons = $right_gene->get_exons();

    foreach my $right_gene_exon (@right_gene_exons) {

        $left_gene->add_mRNA_exon_obj($right_gene_exon);

    }
    
    $left_gene->refine_gene_object();

    return($left_gene);
}

####
sub partition_genes_based_on_gap_overlap {
    my ($all_gaps_aref, $gene_objs_aref) = @_;

    
    my @no_overlap_gap;
    my @with_overlap_gap; # require exon overlap

   
    foreach my $gene (@$gene_objs_aref) {

        my @exon_coords = &get_all_exon_coordinates($gene);

        my $found_overlap = 0;
      overlap_search:
        foreach my $exon_coordset (@exon_coords) {
            my ($exon_lend, $exon_rend) = @$exon_coordset;
            foreach my $gap (@$all_gaps_aref) {
                my ($gap_lend, $gap_rend) = @$gap;

                if ($exon_lend <= $gap_rend && $exon_rend >= $gap_lend) {
                    $found_overlap = 1;
                    last overlap_search;
                }
            }
        }

        if ($found_overlap) {
            push (@with_overlap_gap, $gene);
        }
        else {
            push (@no_overlap_gap, $gene);
        }
        

    }
    
    return(\@no_overlap_gap, \@with_overlap_gap);
    
}

####
sub get_all_exon_coordinates {
    my ($gene_obj) = @_;
    
    my %coords;

    my @exons = $gene_obj->get_exons();

    foreach my $exon (@exons) {

        my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
        $coords{$lend} = $rend;
    }

    my @coord_list;
    
    foreach my $lend (keys %coords) {
        my $rend = $coords{$lend};

        push (@coord_list, [$lend, $rend]);
    }

    return(@coord_list);
}


