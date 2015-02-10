#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 stockholm_alignment  percentColumnOccupied\n\n";

my $stockholm_alignment = $ARGV[0] or die $usage;
my $percent_column_occupied = $ARGV[1] or die $usage;

main: {
	
	
	my @alignments = &parse_alignments($stockholm_alignment);

	## left and right position are zero indexed
	my $left_pos = &get_left_occupied_pos(\@alignments, $percent_column_occupied);

	my $right_pos = &get_right_occupied_pos(\@alignments, $percent_column_occupied);

	my $trimmed_alignment_text = &get_trimmed_alignment_text(\@alignments, $left_pos, $right_pos);

	print $trimmed_alignment_text;
	

	exit(0);
}

####
sub get_trimmed_alignment_text {
	my ($alignments_aref, $left_pos, $right_pos) = @_;

	if ($left_pos > $right_pos) {
		die "Error, left position > right_position!    ($left_pos, $right_pos) ";
	}

	my $alignment_text = "";
	
	foreach my $alignment (@$alignments_aref) {
		my $acc = $alignment->{acc};
		my $alignment_text_aref = $alignment->{alignment_aref};
		
		my $curr_alignment_text = "";
		for (my $i = $left_pos; $i <= $right_pos; $i++) {
			$curr_alignment_text .= $alignment_text_aref->[$i];
		}
		
		$alignment_text .= "$acc $curr_alignment_text\n";

	}

	return ($alignment_text);
}


####
sub get_left_occupied_pos {
	my ($alignments_aref, $min_percent_occupied) = @_;

	my $num_alignments = scalar (@$alignments_aref);
	my $alignment_length = scalar (@{$alignments_aref->[0]->{alignment_aref}});
	
	## examine each alignment position left to right
	for (my $i = 0; $i < $alignment_length; $i++) {
		
		my $num_positions_occupied = 0;
		foreach my $alignment (@$alignments_aref) {
			my $alignment_text_aref = $alignment->{alignment_aref};
			my $char = $alignment_text_aref->[$i];
			if ($char =~ /\w/) {
				$num_positions_occupied++;
			}
		}

		my $percent_occupied = $num_positions_occupied / $num_alignments * 100;

		if ($percent_occupied >= $min_percent_occupied) {
			return ($i);
		}
	}


	## if got here, there's no column position that satisfies the min_percent_occupied criteria
	
	die "Error, no column statisfying the $min_percent_occupied min percent occupied criteria";
	
}

####
sub get_right_occupied_pos {
	my ($alignments_aref, $min_percent_occupied) = @_;

	my $num_alignments = scalar (@$alignments_aref);
	my $alignment_length = scalar (@{$alignments_aref->[0]->{alignment_aref}});
	
	## examine each alignment position right to left
	for (my $i = $alignment_length-1; $i >= 0; $i--) {
		
		my $num_positions_occupied = 0;
		foreach my $alignment (@$alignments_aref) {
			my $alignment_text_aref = $alignment->{alignment_aref};
			my $char = $alignment_text_aref->[$i];
			if ($char =~ /\w/) {
				$num_positions_occupied++;
			}
		}

		my $percent_occupied = $num_positions_occupied / $num_alignments * 100;

		if ($percent_occupied >= $min_percent_occupied) {
			return ($i);
		}
	}


	## if got here, there's no column position that satisfies the min_percent_occupied criteria
	
	die "Error, no column statisfying the $min_percent_occupied min percent occupied criteria";
	
}

	



####
sub parse_alignments {
	my ($alignment_file) = @_;

	my @alignments;

	open (my $fh, $alignment_file) or die "Error, cannot open file $alignment_file";
	while (<$fh>) {
		if (m|//|) {
			# reached end
			last;
		}
		chomp;
		my ($accession, $alignment) = split (/\s+/);
		
		## rid any range indicator
		$accession =~ s|/\d+-\d+||;
		# also remove any : marks, since they interfere with quicktree
		$accession =~ s/:/_/g;
		
		push (@alignments,  { acc => $accession,
							  alignment_aref => [ split (//, $alignment) ],
						  } );
	}
	
	close $fh;

	return (@alignments);
}


							  
		
