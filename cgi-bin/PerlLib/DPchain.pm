#!/usr/bin/env perl

package DPchain;


#####################################################################################################
# all elements in $sorted_eles_aref must be in proper order from left to right
# subroutines:
#    $get_base_score_sref($ele_ref) must return base score for $ele_ref
#    $are_chainable_sref ($eleA_ref, $eleB_ref), where A is before B in list, returns 1 if chainable.
#####################################################################################################

####
sub find_highest_scoring_chain {
	my ($sorted_eles_aref, $get_base_score_sref, $are_chainable_sref) = @_;

	## wrap each ele in a score struct
	my @eles;
	foreach my $ele (@$sorted_eles_aref) {
		my $ele_base_score = &$get_base_score_sref($ele);
		push (@eles, { ele => $ele,
					   base_score => $ele_base_score,
					   sum_score => $ele_base_score,
					   prev => undef,
				   } );
	}

	## do DP chaining:
	
	for (my $i = 1; $i <= $#eles; $i++) {

		my $ele_i = $eles[$i]->{ele};
		my $ele_i_base_score = $eles[$i]->{base_score};
		my $ele_i_sum_score = $eles[$i]->{sum_score};


		for (my $j = $i - 1; $j >= 0; $j--) {
			
			my $ele_j = $eles[$j]->{ele};
			my $ele_j_base_score = $eles[$j]->{base_score};
			my $ele_j_sum_score = $eles[$j]->{sum_score};

			if (&$are_chainable_sref($ele_j, $ele_i) 
				&&
				$ele_j_sum_score + $ele_i_base_score > $ele_i_sum_score) {

				$ele_i_sum_score = $eles[$i]->{sum_score} = $ele_j_sum_score + $ele_i_base_score;
				$eles[$i]->{prev} = $eles[$j];
			}
		}
	}


	## get the highest scoring chain:
	my $best_score = 0;
	my $best_struct = undef;

	foreach my $ele (@eles) {
		if ($ele->{sum_score} > $best_score) {
			$best_score = $ele->{sum_score};
			$best_struct = $ele;
		}
	}

	my @ret_eles;
	while (defined $best_struct) {
		push (@ret_eles, $best_struct->{ele});
		$best_struct = $best_struct->{prev};
	}

	@ret_eles = reverse @ret_eles;
	
	return (@ret_eles);
}


1; #EOM


	
	
