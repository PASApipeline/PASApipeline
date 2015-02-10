#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 pasa.alignment.gff3\n\n";

my $gff3_file = $ARGV[0] or die $usage;


main: {

	my %id_to_features = &parse_alignment_gff3($gff3_file);
	
	foreach my $feature (sort {$a->[0]->{contig} cmp $b->[0]->{contig}} values %id_to_features) {
		
		&write_BED_format($feature);
	}
	
	exit(0);
}

####
sub parse_alignment_gff3 {
	my ($file) = @_;

	my %id_to_features;

	open (my $fh, $file) or die "Error, cannot open file $file";
	while (<$fh>) {
		unless (/\w/) { next; }
		chomp;
		my @x = split(/\t/);
		my $contig = $x[0];
		my $lend = $x[3];
		my $rend = $x[4];
		my $orient = $x[6];
		my $info = $x[8];
		
		my ($id, $target);

		if ($info =~ /ID=([^;]+)/) {
			$id = $1;
		}
		else {
			die "Error, cannot parse ID from $info";
		}

		if ($info =~ /Target=(\S+)/) {
			$target = $1;
		}
		else {
			die "Error, cannot extract target info for $1";
		}

		
		push (@{$id_to_features{$id}}, { lend => $lend,
										 rend => $rend,
										 orient => $orient,
										 target => $target,
										 contig => $contig,

			  });
		
		
	}
	close $fh;
		
	return(%id_to_features);
}


####
sub write_BED_format {
	my ($feature_aref) = @_;

	## compute span
	my ($span_lend, $span_rend);
	## create seg and length lists:
	my @starts;
	my @lengths;
	

	{
		my @coords;
		foreach my $feature (sort {$a->{lend} <=> $b->{lend}} @$feature_aref) {
			
			my $lend = $feature->{lend};
			my $rend = $feature->{rend};
			
			push (@coords, $lend, $rend);
			my $start_pos = $lend - 1;
			my $length = $rend - $lend + 1;
			push (@starts, $start_pos);
			push (@lengths, $length);
			
		}

		@coords = sort {$a<=>$b} @coords;
		

		$span_lend = shift @coords;
		$span_rend = pop @coords;
	}

	my $first_feature = shift @$feature_aref;
	my $contig = $first_feature->{contig};
	my $orient = $first_feature->{orient};
	my $target = $first_feature->{target};

	$span_lend--;
	#$span_rend--; # rend exclusive 
	
	foreach my $start (@starts) {
		$start -= $span_lend;
	}
	
	print join("\t", ($contig, $span_lend, $span_rend, $target, 
					  1, $orient, $span_lend, $span_rend, ".", scalar(@starts),
					  join(",", @lengths),
					  join(",", @starts), 
					  )
		) . "\n";
	
	
	return;
}
	
