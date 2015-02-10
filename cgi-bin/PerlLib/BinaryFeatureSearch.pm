#!/usr/bin/env perl

package BinaryFeatureSearch;


use strict;
use warnings;


sub new {
	my $packagename = shift;

	my $chr_to_features_href = shift;
	
	## Format should be:
	#   chr_to_features_href->{molecule}  = [ a, b, c]
	#  where a = { acc => accession, lend => coord, rend => coord }
		
	my $self = { chr_to_features_href => $chr_to_features_href };

	bless ($self, $packagename);
	
	$self->_sort_features_by_lend();

	return($self);
}

####
sub _sort_features_by_lend {
	my $self = shift;

	my $chr_to_features_href = $self->{chr_to_features_href};
	
	foreach my $molecule (keys %$chr_to_features_href) {
		
		@{$chr_to_features_href->{$molecule}} = sort {$a->{lend}<=>$b->{lend}} @{$chr_to_features_href->{$molecule}};
	}

	return;
}



####
sub find_overlapping_feature {
	my $self = shift;
	my ($molecule, $pos_lend, $pos_rend) = @_;
	
	my $feature_list_aref = $self->{chr_to_features_href}->{$molecule};
	unless (ref $feature_list_aref) {
		print STDERR "No features for $molecule\n";
		return;
	}

	## perform binary search
	
	my $i_left = 0; 
	my $i_right = $#$feature_list_aref;

	
	
	while (1) {
		
		if ($i_right < $i_left) {
			last;
		}

		my $mid_index = int(($i_left+$i_right)/2 + 0.5);
		
		my $feature = $feature_list_aref->[$mid_index];
		
		my ($acc, $lend, $rend) = ($feature->{acc},
								   $feature->{lend},
								   $feature->{rend});
		
		# print "Comparing [$i_left, $i_right, $mid_index] to $acc\t$lend\t$rend\n";

		if ($pos_rend < $lend) {
			# binary search left
			$i_right = $mid_index - 1;
			next;
		}
		elsif ($pos_lend > $rend) {
			$i_left = $mid_index + 1;
			next;
		}
		elsif ($lend <= $pos_rend && $rend >= $pos_lend) {
			return($feature);
		}
		else {
			# not sure what to do:
			die "Error, doesn't overlap and isn't on either side";
		}
	}
	return;
}

1; #EOM
