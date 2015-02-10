#!/usr/bin/env perl

use strict;
use warnings;
use lib ($ENV{EUK_MODULES});
require "overlapping_nucs.ph";


my $MAX_INTRON_LEN = 1000; #should make this a variable in the future.

my $usage = "usage: $0 tab_delim_summary_file (m2fmt|btab)\n\n";

my $input_file = $ARGV[0] or die $usage;
my $fmt = $ARGV[1] or die $usage;

unless ($fmt =~ /m2fmt|btab/) { 
	die $usage;
}

my $MAX_PERCENT_HIT_OVERLAP = 20;

main: {

	my @hits;
	
	if ($fmt =~ /m2fmt/) {
		@hits = &parse_hits_from_m2file($input_file);
	}
	else {
		@hits = &parse_hits_from_btab_file($input_file);
	}

	if (@hits) {
		
		my @genome_seq_partitioned_hits = &partition_hits_by_contig(@hits);

		foreach my $hits_aref (@genome_seq_partitioned_hits) {
			
			my @chains = &chain_hits (@$hits_aref);
		
			&print_chains(@chains);
		}
	}
	
	exit(0);
}

####
sub partition_hits_by_contig {
	my @hits = @_;


	my %contig_to_hit_list;

	foreach my $hit (@hits) {
		my $contig = $hit->{sid};
		push (@{$contig_to_hit_list{$contig}}, $hit);
	}


	return (values %contig_to_hit_list);
}

####
sub parse_hits_from_btab_file {
	my ($btab_file) = @_;
	
	my @hits;
	
	open (my $fh, $btab_file) or die $!;
	while (<$fh>) {
		my $line = $_;
		chomp;
		my @x = split (/\t/);
		
		my ($qid, $sid, $query_end5, $query_end3, $search_end5, $search_end3) = ($x[0], $x[5], $x[6], $x[7], $x[8], $x[9]);
		
		my $per_id = $x[10];
		

		my $orient = ($search_end5 < $search_end3) ? '+' : '-';
		
		my $query_genome_acc = "$qid;;$sid";
		
		my $hit_struct = { qid => $qid,
						   sid => $sid,
						   
						   qstart => $query_end5,
						   qend => $query_end3,
						   sstart => $search_end5,
						   send => $search_end3,
						   
						   sorient => $orient,
						   
						   line => $line,
						   						   
						   query_genome_acc => $query_genome_acc,
						   
						   per_id => $per_id,
						   
					   };

		push (@hits, $hit_struct);
	}


	return (@hits);
}


####
sub parse_hits_from_m2file {
	my ($m2file) = @_;

	my @hits;

	open (my $fh, $m2file) or die "Error, cannot open file $m2file";
	while (<$fh>) {
		my $line = $_;
		chomp;
		my ($qid,$sid,$E, $N, $Sprime, $S, $alignlen, $nident, $npos,
			$nmism, $pcident, $pcpos, $qgaps, $qgaplen, $sgaps, $sgaplen,
			$qframe, $qstart, $qend, $sframe, $sstart, $send) = split (/\t/);
		

		my $query_genome_acc = "$qid;;$sid";
		
		my $sorient = ($sframe =~ /\+/) ? '+' : '-';
		
		my $hit_struct = { qid => $qid,
						   sid => $sid,
						   
						   qstart => $qstart,
						   qend => $qend,
						   sstart => $sstart,
						   send => $send,

						   sorient => $sorient,
	
						   line => $line,

						   query_genome_acc => $query_genome_acc, # used to simplify query and genome accession groupings.

						   per_id => $pcident,
						   
					   };
		
		
		push (@hits, $hit_struct);
	}
	close $fh;
	
	unless (@hits) {
		die "No hits parsed from file $m2file";
	}
	
	return (@hits);
}


####
sub chain_hits {
	my @hits = @_;

	## sort by genome coordinate
	@hits = sort {
		
		$a->{query_genome_acc} cmp $b->{query_genome_acc}
	  	     ||
			$a->{sstart}<=>$b->{sstart}

	} @hits;

	
	my @chains = ( [ $hits[0] ] );

	for (my $i = 1; $i <= $#hits; $i++) {

		my $last_chain_aref = $chains[$#chains];


		my $prev_hit = $last_chain_aref->[ $#{$last_chain_aref} ];
		my ($prev_hit_qLEND, $prev_hit_qREND) = sort {$a<=>$b} ($prev_hit->{qstart}, $prev_hit->{qend});
		my ($prev_hit_sLEND, $prev_hit_sREND) = sort {$a<=>$b} ($prev_hit->{sstart}, $prev_hit->{send});

		my $curr_hit = $hits[$i];
		my ($curr_hit_qLEND, $curr_hit_qREND) = sort {$a<=>$b} ($curr_hit->{qstart}, $curr_hit->{qend});
		my ($curr_hit_sLEND, $curr_hit_sREND) = sort {$a<=>$b} ($curr_hit->{sstart}, $curr_hit->{send});
		
		if ( $prev_hit->{query_genome_acc} eq $curr_hit->{query_genome_acc}
			 &&
			 $prev_hit->{sorient} eq $curr_hit->{sorient}
			 &&
			 (! &overlap_too_much([$prev_hit_qLEND, $prev_hit_qREND], [$curr_hit_qLEND, $curr_hit_qREND]))
			 &&
			 
			 # tentative intron length check
			 (! &too_far_apart([$prev_hit_sLEND, $prev_hit_sREND], [$curr_hit_sLEND, $curr_hit_sREND])) ## or too far apart!
			 
			 &&
			 (
			 ($curr_hit->{sorient} eq '+' && $curr_hit_qREND > $prev_hit_qREND) 
			                              ||
			  ($curr_hit->{sorient} eq '-' && $curr_hit_qREND < $prev_hit_qREND) 
			  )
			 ) 
		{
			
			## append to current chain:
			push (@$last_chain_aref, $curr_hit);
		}
		else {
			# create a new chain:
			push (@chains, [ $curr_hit ] );
		}
		

	}
		
	return (@chains);
}
		
####
sub overlap_too_much {
	my ($coordset_A_aref, $coordset_B_aref) = @_;

	unless (&coordsets_overlap($coordset_A_aref, $coordset_B_aref)) {
		return (0);
	}

	my ($A1, $A2) = sort {$a<=>$b} @$coordset_A_aref;
	my ($B1, $B2) = sort {$a<=>$b} @$coordset_B_aref;

	my $A_len = abs ($A2 - $A1) + 1;
	my $B_len = abs ($B2 - $B1) + 1;

	my $overlap = &nucs_in_common($A1, $A2, $B1, $B2);
	
	my $percent_overlap_A = $overlap / $A_len * 100;
	my $percent_overlap_B = $overlap / $B_len * 100;

	if ($percent_overlap_A > $MAX_PERCENT_HIT_OVERLAP || $percent_overlap_B > $MAX_PERCENT_HIT_OVERLAP) {
		## coordinate pairs exceed allowed overlap
		return (1);
	}
	else {
		# all good!
		return (0);
	}
}


####
sub too_far_apart {
	my ($coordset_A_aref, $coordset_B_aref) = @_;

	if (&coordsets_overlap($coordset_A_aref, $coordset_B_aref)) {
		return (0);
	}

	my ($A1, $A2) = sort {$a<=>$b} @$coordset_A_aref;
	my ($B1, $B2) = sort {$a<=>$b} @$coordset_B_aref;

	if (abs($B1 - $A2) -1 > $MAX_INTRON_LEN) {
		return (1); # Too far apart!
	}
	else {
		return (0); # within range.
	}
}



	
####
sub print_chains {
	my @chains = @_;

	foreach my $chain (@chains) {
		
		my $chain_text = "";

		my $sum_score = 0;

		my $orient;
		my @genome_coords;
		my @prot_coords;
		my ($qid, $sid);

		
		foreach my $hit (@$chain) {
			
			$chain_text .= $hit->{line};
			
			my $per_id = $hit->{per_id};
			my ($qstart, $qend) = ($hit->{qstart}, $hit->{qend});
			my $hit_len = abs ($qend - $qstart) + 1;
			$sum_score += $hit_len * $per_id / 100;
			
			$orient = $hit->{sorient};
			$qid = $hit->{qid};
			$sid = $hit->{sid};
			
			push (@genome_coords, $hit->{sstart}, $hit->{send});
			push (@prot_coords, $qstart, $qend);

		}
	
		@genome_coords = sort {$a<=>$b} @genome_coords;
		@prot_coords = sort {$a<=>$b} @prot_coords;
		
		my $genome_lend = shift @genome_coords;
		my $genome_rend = pop @genome_coords;
		
		my $prot_lend = shift @prot_coords;
		my $prot_rend = pop @prot_coords;

		my $range_text = "$qid\t$prot_lend-$prot_rend\t$sid\t$genome_lend-$genome_rend\t$orient";

		print "\n#Chain\t$range_text\t" . sprintf("%.1f", $sum_score) . "\n" . $chain_text;
		
	}
	
	return;
}


		
