#!/usr/bin/env perl


use strict;
use warnings;
use Data::Dumper;

my $usage = "usage: $0 gene_annotations.gff3\n\n";

my $gff3_file = $ARGV[0] or die $usage;


my %feature_ID_to_data_struct;
my %parent_to_children;
my %child_to_parent;

my %gene_component_hierarchy = ( 'mRNA' => 'gene',
								 'exon' => 'mRNA',
								 'CDS' => 'mRNA',
								 );

my $has_ERROR_flag = 0;

main: {

  &parse_GFF3_file($gff3_file);

  &check_child_parent_consistency();

  &check_features_have_parents();

  &ensure_CDS_and_exon_encapsulation();

  print "\n\nDone validating record. Warnings would have been reported if potential problems were detected.\n\n";
  
  exit($has_ERROR_flag);
  
}


####
sub parse_GFF3_file {
    my ($gff3_file) = @_;
    
    open (my $fh, $gff3_file) or die "Error, cannot open file $gff3_file";
    while (<$fh>) {
        chomp;
        if (/^\#/) { next; } # comment line
        unless (/\w/) { next;}
        
        my @x = split (/\t/);
        
        my ($contig_id, $feat_type, $lend, $rend, $orient, $feature_info) = ($x[0], $x[2], $x[3], $x[4], $x[6], $x[8]);
        
        ## the only features we care about are gene, mRNA, exon, and CDS
        unless ($feat_type =~ /^(gene|mRNA|exon|CDS)$/) { next; }
        
        $feature_info =~ /ID\s*=\s*([^;]+)/ or die "Fatal Error: cannot parse ID from entry\n$_";
        
        my $feature_ID = $1;
        $feature_ID =~ s/\s+$//;
        
        my $parent_ID = undef;
        
        if ($feature_info =~ /Parent\s*=\s*([^;]+)/) {
            $parent_ID = $1;
            $parent_ID =~ s/\s+$//;
        }
        
        my $feature_struct = { feat_type => $feat_type,
                               lend => $lend,
                               rend => $rend,
                               orient => $orient,
                               contig => $contig_id,
                               parent_ID => $parent_ID,
                               feature_ID => $feature_ID
                               };	
        
        if (exists ($feature_ID_to_data_struct{$feature_ID})) {
            ## allow it for now, but make sure the data are the same!
            my $stored_feature = $feature_ID_to_data_struct{$feature_ID};
            if (! identical_features($stored_feature, $feature_struct)) {
                print STDERR "Error, feature: $feature_ID is described multiple times with different data values:\n" . Dumper($feature_struct) . Dumper ($stored_feature);
                $has_ERROR_flag = 1;
                next;
            }
            else {
                ## no new info here. Already stored it, so go on to next entry
                next;
            }
            
        }
        else {
            # store it:
            $feature_ID_to_data_struct{$feature_ID} = $feature_struct;
        }
        
        if ($parent_ID) {
            push (@{$parent_to_children{$parent_ID}}, $feature_ID);
            $child_to_parent{$feature_ID} = $parent_ID;
        }
        
        
    }
    
    return;
}

####
sub check_child_parent_consistency {
    
    ## check:
    ## -same contig ID
    ## -same orientation
    ## -child coordinates are encapsulated by parents
    ## -valid parent/child relationship according to component hierarchy
    
    
    foreach my $parent_ID (keys %parent_to_children) {
        
        my $parent_struct = $feature_ID_to_data_struct{$parent_ID};
        
        unless (ref $parent_struct) {
            die "Fatal Error, cannot locate data entry for ID: [$parent_ID]";
        }
        my ($feat_type, $lend, $rend, $orient, $contig) = ($parent_struct->{feat_type},
                                                           $parent_struct->{lend},
                                                           $parent_struct->{rend},
                                                           $parent_struct->{orient},
                                                           $parent_struct->{contig});
        
        
        
        my @children_IDs = @{$parent_to_children{$parent_ID}};
        
        foreach my $child (@children_IDs) {
            
            my $child_struct = $feature_ID_to_data_struct{$child} or die "Fatal error, cannot locate data entry for ID: $child";
            my ($child_feat_type, $child_lend, $child_rend, $child_orient, $child_contig) = ($child_struct->{feat_type},
                                                                                             $child_struct->{lend},
                                                                                             $child_struct->{rend},
                                                                                             $child_struct->{orient},
                                                                                             $child_struct->{contig});
            
            
            ## check contig ID
            if ($contig ne $child_contig) {
                print "Error, parent $parent_ID ($feat_type, $contig) and child $child ($child_feat_type, $child_contig) are located on different contigs\n";
                $has_ERROR_flag = 1;
            }
            
            ## check orientation:
            if ($orient ne $child_orient) {
                print "Error, parent $parent_ID ($feat_type, $orient) and child $child ($child_feat_type, $child_orient) have conflicting orientations. \n";
                $has_ERROR_flag = 1;
            }
            
            ## check coordinate encapsulation:
            if (!  ($child_lend >= $lend && $child_rend <= $rend) ) {
                
                print "Error, parent $parent_ID ($feat_type, $lend-$rend) does not encapsulate coords of child $child ($child_feat_type, $child_lend-$child_rend) \n";
                $has_ERROR_flag = 1;
            }
            
            ## check for hierarchy compatibility
            unless ($gene_component_hierarchy{$child_feat_type} eq $feat_type) {
                print "Error, parent $parent_ID ($feat_type) cannot have a child $child of type $child_feat_type\n";
                $has_ERROR_flag = 1;
            }
            
        }
        
		
    }
    
    return;
}


####
sub identical_features {
    my ($featureA, $featureB) = @_;
    
    if ($featureA->{feat_type} eq 'CDS' && $featureB->{feat_type} eq 'CDS' &&
        $featureA->{parent_ID} eq $featureB->{parent_ID}) {
        
        ## same CDS identifier, that's allowed
        return(1);
    }
    
    
    foreach my $feat_key (keys %$featureA) {
        if (defined $featureA->{$feat_key} && defined $featureB->{$feat_key} 
            && $featureA->{$feat_key} ne $featureB->{$feat_key}) {
            return (0);
        }
        
    }
    if (defined($featureA->{parent_ID}) xor defined ($featureB->{parent_ID}) ) {
        return (0);
    }
    
    ## if got this far, then all populated data fields should be identical
    return (1); 
}

####
sub ensure_CDS_and_exon_encapsulation {
    
    ## examine every mRNA and the linked exons and CDSs.
    ## Every mRNA *must* have at least one CDS record, as required by PASA.
    ## Also, each CDS should be encapsulated by a single exon
    
    my @mRNAs;
    foreach my $data_struct (values %feature_ID_to_data_struct) {
        if ($data_struct->{feat_type} eq 'mRNA') {
            my $mRNA_feature_ID = $data_struct->{feature_ID};
            my @children_IDs = @{$parent_to_children{$mRNA_feature_ID}};
            
            my @exons;
            my @CDSs;
            foreach my $child_ID (@children_IDs) {
                my $child_data_struct = $feature_ID_to_data_struct{$child_ID};
                if ($child_data_struct->{feat_type} eq 'exon') {
                    push (@exons, $child_data_struct);
                }
                elsif ($child_data_struct->{feat_type} eq 'CDS') {
                    push (@CDSs, $child_data_struct);
                }
                else {
                    die "Error, found a feature as a child to an mRNA with an unsupported child feature type: " . Dumper ($child_data_struct);
                }
            }
            
            unless (@CDSs) {
                print "Warning: mRNA $mRNA_feature_ID lacks a CDS record - OK if a ncRNA gene.\n";
            }
            &correlate_CDSs_with_exons(\@exons, \@CDSs);
            
        }
    }
    
    return;
}

####
sub correlate_CDSs_with_exons {
    my ($exons_aref, $CDSs_aref) = @_;
    
    my %exon_used;
    
    foreach my $CDS (@$CDSs_aref) {
        
        my $cds_ID = $CDS->{feature_ID};
        my ($cds_lend, $cds_rend) = ($CDS->{lend}, $CDS->{rend});
        
        ## map to an exon:
        my $found_exon = 0;
        foreach my $exon (@$exons_aref) {
            my $exon_ID = $exon->{feature_ID};
            my ($exon_lend, $exon_rend) = ($exon->{lend}, $exon->{rend});
            
            
            if ($cds_lend >= $exon_lend && $cds_rend <= $exon_rend) {
                # encapsulated:
                if (my $other_ID = $exon_used{$exon_ID}) {
                    print "ERROR, CDS $cds_ID ($cds_lend-$cds_rend) maps to exon $exon_ID ($exon_lend-$exon_rend), but this exon already encodes a different CDS record $other_ID\n";
                    $has_ERROR_flag = 1;
                }
                
                $found_exon=1;
                $exon_used{$exon_ID} = $cds_ID;
                last;
            }
        }
        
        if (! $found_exon) {
            print "ERROR, CDS $cds_ID does not fully map within an exon record.\n";
            $has_ERROR_flag = 1;
        }
        
    }
    
    return;
}

####
sub check_features_have_parents {
    
    foreach my $feature_ID (keys %feature_ID_to_data_struct) {
        
        my $struct = $feature_ID_to_data_struct{$feature_ID};
        my $feat_type = $struct->{feat_type};
        
        if ($feat_type ne 'gene') {
            # must have a parent!
            
            my $parent_feat_ID = $child_to_parent{$feature_ID};
            unless (defined $parent_feat_ID) {
                print "ERROR, feature $feature_ID ($feat_type) lacks a parent feature\n";
                $has_ERROR_flag = 1;
                next;
            }
            
            my $parent_struct = $feature_ID_to_data_struct{$parent_feat_ID};
            my $parent_feat_type = $parent_struct->{feat_type};
            unless ($gene_component_hierarchy{$feat_type} eq $parent_feat_type) {
                print "ERROR, feature $feature_ID ($feat_type) has a parent $parent_feat_ID ($parent_feat_type) and is not allowed!\n";
                $has_ERROR_flag = 1;
            }
        }
        
        
    }
    
    return;
}

