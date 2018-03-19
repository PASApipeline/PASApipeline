## This package defines a namespace only.  DO NOT INSTANTIATE

## Does NOT operate on isoforms.  Each isoform must be submitted individually.

our $DEBUG;
our $SEE;
package Gene_obj_merger;
use strict;
use Gene_obj;


my $temp_id = 0; #used to assign temp_ids to incoming gene components.

sub merge_gene_models {
    my ($gene_obj1, $gene_obj2) = @_;
    print "Merging gene models\n";
    ## assign temporary feat_names to uniquely identify the gene components (EXON, CDS).
    &assign_exon_temp_ids($gene_obj1);
    &assign_exon_temp_ids($gene_obj2);

    my $cloned_gene = $gene_obj1->clone_gene();; #new gene to be the result of the merge.
    $cloned_gene->{alt_locus} = $gene_obj2->{locus};
    
    ## build cloned_gene gene structure based on merge.
    $cloned_gene->erase_gene_structure();
    $cloned_gene->delete_isoforms();
    my %merged; # track feat_names that are merged.
    my @gene_obj1_exons = $gene_obj1->get_exons();
    my @gene_obj2_exons = $gene_obj2->get_exons();
    foreach my $exon1 (@gene_obj1_exons) {
	foreach my $exon2 (@gene_obj2_exons) {
	    if (&are_diradj_or_overlap($exon1, $exon2)) {
		print "Exons overlap and are being merged now\n" if $DEBUG||$SEE;
		print "merging:\n " . $exon1->toString() . $exon2->toString() if $DEBUG||$SEE; 
		$cloned_gene->add_mRNA_exon_obj (&merge_exons($exon1, $exon2));
		$merged{$exon1->{temp_id}} = 1;
		$merged{$exon2->{temp_id}} = 1;
	    }
	}
    }
    ## add unmerged exons to gene model
    foreach my $exon (@gene_obj1_exons, @gene_obj2_exons) {
	unless ($merged{$exon->{temp_id}}) {
	    print "Adding unmerged exon to gene\n" if $DEBUG||$SEE;
	    print "adding:\n" . $exon->toString() if $DEBUG||$SEE;
	    $cloned_gene->add_mRNA_exon_obj($exon->clone_exon());
	}
    }
    $cloned_gene->refine_gene_object();
    return ($cloned_gene);
}


####
sub are_diradj_or_overlap {
    #determines if the two coordinate sets are directly adjacent or overlap each other.
    my ($exon1, $exon2) = @_;
    my @x = $exon1->get_mRNA_exon_end5_end3();
    ## put all coord pairs in a positive direction.
    my @coordpair1 = ($x[0]<$x[1]) ? @x : (reverse @x);
    # offset coord pairs so can look for overlap instead of adjacent.
    $coordpair1[0]--;
    $coordpair1[1]++;
    ## now do the same with the other exon
    my @y = $exon2->get_mRNA_exon_end5_end3();
    my @coordpair2 = ($y[0]<$y[1]) ? @y : (reverse @y);
    $coordpair2[0]--;
    $coordpair2[1]++;
    # test for overlap
    if ($coordpair1[0] <= $coordpair2[1] && $coordpair1[1] >= $coordpair2[0]) {
	return (1); #true
    } else {
	return (0); # false
    }
}



#### 
sub merge_exons {
    my ($exon1, $exon2) = @_;
    my @exon_coords = (sort {$a<=>$b} ($exon1->get_mRNA_exon_end5_end3(), $exon2->get_mRNA_exon_end5_end3()));
    my @cds_coords = (sort {$a<=>$b} ($exon1->get_CDS_end5_end3(), $exon2->get_CDS_end5_end3()));
    my ($exon_end5, $exon_end3, $cds_end5, $cds_end3); #merged exon values.
    my $orientation = $exon1->get_orientation();
    # assume forward orientation
    ## process rna exon info
    ($exon_end5, $exon_end3) = ($exon_coords[0], $exon_coords[$#exon_coords]);
    if ($orientation eq '-') {
	($exon_end5, $exon_end3) = ($exon_end3, $exon_end5);
    }
    my $new_exon = new mRNA_exon_obj ($exon_end5, $exon_end3);
    ## process cds-exon info
    if (@cds_coords) {
	($cds_end5, $cds_end3) = ($cds_coords[0], $cds_coords[$#cds_coords]);
	if ($orientation eq '-') {
	    ($cds_end5, $cds_end3) = ($cds_end3, $cds_end5);
	}
	$new_exon->add_CDS_exon_obj($cds_end5, $cds_end3);
    }
    return ($new_exon);
}

####
sub assign_exon_temp_ids {
    my $gene_obj = shift;
    my @exons = $gene_obj->get_exons();
    foreach my $exon (@exons) {
	$temp_id++;
	$exon->{temp_id} = $temp_id;
	my $cds = $exon->get_CDS_obj();
	if (ref $cds) {
	    $temp_id++;
	    $cds->{temp_id} = $temp_id;
	}
    }
}

    
1; #return success.
