#!/usr/local/bin/perl

## This module is used to compare to Gene_objs.  Everything here is 'static'.  No intentional object members.
## the static fields can be reset by calling the reset method.
## The major assumption here is that the gene_obj's are on the same sequence.  Only coordinates are compared, not seq info.
## used to compare Ceres cDNA alignment-based gene models to current Arabidopsis annotations.

package main;
our $SEE;

=head1 NAME

Gene_obj_comparator.pm

=cut

=head1 DESCRIPTION

This module is used to compare two gene objects (Gene_obj.pm) to see if they are identical, have the same CDS features, and/or same mRNA features.  This package does not need to be instantiated.  The methods are exported so they do not have to be qualified by the package name.

First, use the compare_genes() method to perform the comparison, then use the other public methods to capture the results of the comparison.

=cut




package Gene_obj_comparator;
require Exporter;
use strict;
our @ISA = qw(Exporter);
our @EXPORT = qw(compare_genes are_CDS_same are_mRNA_same reset are_identical analyze_cds_frames);


## file scoped lexicals.  access via access methods.  Reset to false via reset method.
my $CDS_same = 0; #initialize to false.
my $mRNA_same = 0; #initialize to false.

=item compare_gene_structures()

=over 4

B<Description:> This method performs the comparison between two Gene_obj objects.  This method should be performed first, followed by the methods below to ascertain similiarities/differences.

B<Parameters:> $gene_obj1, $gene_obj2

$gene_obj1 and $gene_obj2 are of type Gene_obj specified in Gene_obj.pm

B<Returns:> [1|0]

This method returns the truth/false of an identity test.  Identity tests are performed only based on coordinates.

** No alternative splicing support **  Only comparisons are done at the primary gene level.  If alt isoforms exist, they should be compared separately.

=back

=cut




## module is now extended for functional annotation comparisons as well.
sub compare_gene_structures {
    return (&compare_genes(@_));
}


## old method, left in for backwards compatibility.  Called now by compare_gene_structures().
sub compare_genes {
    my ($gene1, $gene2) = @_;
    unless (ref ($gene1) && ref($gene2)) {
	print STDERR "I require two gene objects to compare!\n";
	return (undef());
    }

    &reset();

    ## Gather coord info
    my (%CDS_coords1, %CDS_coords2, %mRNA_coords1, %mRNA_coords2);
    &populate_gene_data ($gene1, \%mRNA_coords1, \%CDS_coords1);
    &populate_gene_data($gene2, \%mRNA_coords2, \%CDS_coords2);
    
    ## perform coord comparison
    $mRNA_same = &compare_coordsets(\%mRNA_coords1, \%mRNA_coords2);
    $CDS_same = &compare_coordsets (\%CDS_coords1, \%CDS_coords2);
    return ($mRNA_same && $CDS_same);
}
    

=item are_CDS_same()

=over 4

B<Description:> Determines whether two gene objects have identical CDS structures.

CDS refers to the protein coding segments of individual exons.

B<Parameters:> none

B<Returns:> [1|0]

1 (true) is returned if the CDS are the same.

0 (false) is returned if there are any differences in CDS structures.

=back

=cut




sub are_CDS_same {
    return ($CDS_same);
}


=item are_mRNA_same()

=over 4

B<Description:> Determines if the mRNA exon structure is identical between two gene objects.

The mRNA exons are the intron-less segments of the full-length transcript, including the CDS features and untranslated regions (UTRs).

B<Parameters:> none.

B<Returns:> [1|0]

1=true, 0=false.

=back

=cut


sub are_mRNA_same {
    return ($mRNA_same);
}



=item are_identical()

=over 4

B<Description:> Provides the result of an identity test between two objects.  Identity tests are performed only at the coordinate level.  An identical set of gene objects have identical mRNA exons and identical CDS structure coordinates.

B<Parameters:> none.

B<Returns:> [0|1]

1=true, 0=false.

=back

=cut




sub are_identical {
    return ($CDS_same && $mRNA_same);
}


sub reset {
    $CDS_same = 0;
    $mRNA_same = 0;
}


####
sub compare_functional_annotations {
    my ($A_gene, $B_gene) = @_;
    unless (ref $A_gene && ref $B_gene) {
	die "Must provide two Gene_obj references.\n";
    }
    
    my @updates;
    
    ## check for different com_name annotation.
    my $a_com_name = $A_gene->get_product_names();
    my $b_com_name = $B_gene->get_product_names();
    if ($a_com_name ne $b_com_name) {
	push (@updates, "product_name_change");
    }
    
    ## check for different pseudogene annotation.
    if ($A_gene->{is_pseudogene} != $B_gene->{is_pseudogene}) {
	push (@updates, "pseudogene_change");
    }
    
    ## check for gene symbol change:
    if ($A_gene->get_gene_symbols() ne $B_gene->get_gene_symbols()) {
	push (@updates, "gene_symbol_change");
    }
    
    ## check for ec number change
    if ($A_gene->get_ec_numbers() ne $B_gene->get_ec_numbers() ) {
	push (@updates, "ec_number_change");
    }
    
    ## check for gene name change
    if ($A_gene->get_gene_names() ne $B_gene->get_gene_names()) {
	push (@updates, "gene_name_change");
    }
    
    ## Check for Gene Ontology updates:
    my $same_GO = &same_GO_assignments ($A_gene, $B_gene);
    unless ($same_GO) {
	push (@updates, "GeneOntology_change");
    }

    return (@updates);
    
}


######################
# private methods
 
sub populate_gene_data {
    my ($gene_obj, $mRNA_coords_ref, $CDS_coords_ref) = @_;
    my @exons = $gene_obj->get_exons();
    foreach my $exon (@exons) {
	my ($mRNA_end5, $mRNA_end3) = $exon->get_mRNA_exon_end5_end3();
	$mRNA_coords_ref->{$mRNA_end5} = $mRNA_end3;
	my @CDS_coords = $exon->get_CDS_end5_end3();
	if (@CDS_coords) {
	    my ($cds_end5, $cds_end3) = @CDS_coords;
	    $CDS_coords_ref->{$cds_end5} = $cds_end3;
	}
    }
}


sub compare_coordsets {
    my ($cset1_ref, $cset2_ref) = @_;
    my $same = 1;
    foreach my $end5 (keys %$cset1_ref, keys %$cset2_ref) {
	if ($cset1_ref->{$end5} != $cset2_ref->{$end5}) {
	    $same = 0;
	    last;
	}
    }
    return ($same);
}



####
sub same_GO_assignments {
    my ($A_gene, $B_gene) = @_;
    
    my @A_gene_GO_assignments = $A_gene->get_gene_ontology_objs();
    my $A_go_text = "";
    my @A_go_texts;
    foreach my $go_assignment (@A_gene_GO_assignments) {
	push (@A_go_texts, $go_assignment->toString());
    }
    $A_go_text = join ("\n", @A_go_texts);
    
    my @B_gene_GO_assignments = $B_gene->get_gene_ontology_objs();
    my $B_go_text = "";
    my @B_go_texts;
    foreach my $go_assignment (@B_gene_GO_assignments) {
	push (@B_go_texts, $go_assignment->toString());
    }
    $B_go_text = join ("\n", @B_go_texts);
    
    if ($A_go_text eq $B_go_text) {
	return (1);
    } else {
	return (0);
    }
}




=item analyze_cds_frames()

=over 4

B<Description:> Given two gene objects (geneA, geneB), the coding region and introns of 
                geneA are compared to the coding region and introns of geneB

B<Parameters:> (geneA, geneB)

B<Returns:> compare_struct_href


    compare_struct_href = {
        
        coding_coding_same => 0|1,  CDS regions overlap and are in the same frame
        
        coding_coding_diff => 0|1,  CDS regions overlap and are in different frames
        
        coding_intron => 0|1,       geneA coding region overlaps geneB intron
        
        intron_coding => 0|1,       geneA intron overlaps geneB coding
                       
        intron_intron => 0|1,       geneA intron overlaps geneB intron
        
    };



=back

=cut



####
sub analyze_cds_frames {
    my ($gene_obj_A, $gene_obj_B) = @_;
    
    my @coords = ($gene_obj_A->get_gene_span(), $gene_obj_B->get_gene_span());
    
    @coords = sort {$a<=>$b} @coords;

    my $min_lend = shift @coords;
    my $max_rend = pop @coords;

    my $length = $max_rend - $min_lend + 1;

    ## store frame for each bp.
    
    #  extragenic:                (-)
    #  intron:                    (.)
    #  cds:                     ([123])

    my @geneA_frames;
    my @geneB_frames;
    
    &_build_structure_vector (\@geneA_frames, $length, $min_lend, $gene_obj_A);
    &_build_structure_vector (\@geneB_frames, $length, $min_lend, $gene_obj_B);
        
    ## just looking witin the model span of the template, ignoring utrs.
    
    my $coding_coding_same = 0;
    my $coding_coding_diff = 0;
    my $coding_intron = 0;
    my $intron_coding = 0;
    my $intron_intron = 0;
    
    for (my $i = 0; $i <= $length; $i++) {
        
        my $geneA_class = $geneA_frames[$i];
        my $geneB_class = $geneB_frames[$i];
        

        if ($geneA_class ne '-' && $geneB_class ne '-') {
            if ($geneA_class eq '.' || $geneB_class eq '.') {
                
                ## at least one is in the intron
                if ($geneA_class eq '.' && $geneB_class eq '.') {
                    $intron_intron = 1;
                }
                elsif ($geneA_class ne '.') {
                    $coding_intron = 1;
                }
                elsif ($geneB_class ne '.') {
                    $intron_coding = 1;
                }
            }
            else {
                ## neither in intron, so must be coding/coding comparison
                if ($geneA_class eq $geneB_class) {
                    $coding_coding_same = 1;
                }
                else {
                    # must be in different frames
                    $coding_coding_diff = 1;
                }
            }
        }
    }
    
    my $ret_struct = { coding_coding_same => $coding_coding_same,
                       coding_coding_diff => $coding_coding_diff,
                       
                       coding_intron => $coding_intron,
                       intron_coding => $intron_coding,
                       
                       intron_intron => $intron_intron,
                   };


    return ($ret_struct);
}



sub _build_structure_vector {
    my ($structure_aref, $vector_length, $lend_adjust, $gene_obj) = @_;
    
    #  extragenic:                (-)
    #  intron:                    (.)
    #  cds:                     ([123])
    

    # init all as extragenic
    for (my $i = 0; $i <= $vector_length; $i++) {
        $structure_aref->[$i] = "-";
    }
    
    
    ## build coding vector portion
    foreach my $exon ($gene_obj->get_exons()) {
        my $cds = $exon->get_CDS_obj();
        if (ref $cds) {
            my ($end5, $end3) = $cds->get_coords();
                        
            $end5 -= $lend_adjust;
            $end3 -= $lend_adjust;
            
            my $phase = $cds->{phase};
            unless (defined $phase) {
                print STDERR "Warning, phase is not set for cds." . $gene_obj->toString() . "\n\t** defaulting to zero.";
                $phase = 0;
            }
            
            my $count = 1 + $phase;
            if ($gene_obj->get_orientation() eq '+') {
                for (my $i = $end5; $i <= $end3; $i++) {
                    $structure_aref->[$i] = ($count - 1) % 3;
                    $count++;
                }
            } 
            else {
                # minus strand
                for (my $i = $end5; $i >= $end3; $i--) {
                    $structure_aref->[$i] = ($count - 1) % 3;
                    $count++;
                }
            }
        }
    }
    
    
    ## build intron vector portion
    foreach my $intron_coordset ($gene_obj->get_intron_coordinates()) {
        my ($intron_lend, $intron_rend) = sort {$a<=>$b} @$intron_coordset;
        
        $intron_lend -= $lend_adjust;
        $intron_rend -= $lend_adjust;

        for (my $i = $intron_lend; $i <= $intron_rend; $i++) {
            $structure_aref->[$i] = '.';
        }
    }


}



1; #end of module
