#!/usr/bin/env perl

use strict;
use warnings;
use lib ($ENV{EUK_MODULES});
require "overlapping_nucs.ph";

my $usage = "usage: $0 m2formatTblastNresultFile\n\n";

my $m2file = $ARGV[0] or die $usage;


my $MAX_PERCENT_HIT_OVERLAP = 20;

main: {

	my @hits = &parse_hits_from_m2file($m2file);
	
	my @chains = &chain_hits (@hits);
	
	&print_chains(@chains);


	exit(0);
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
						   sframe => $sframe,
						   qstart => $qstart,
						   qend => $qend,
						   sstart => $sstart,
						   send => $send,

						   sorient => $sorient,
	
						   line => $line,

						   query_genome_acc => $query_genome_acc, # used to simplify query and genome accession groupings.
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
		

		my $curr_hit = $hits[$i];
		my ($curr_hit_qLEND, $curr_hit_qREND) = sort {$a<=>$b} ($curr_hit->{qstart}, $curr_hit->{qend});
				
		if ( $prev_hit->{query_genome_acc} eq $curr_hit->{query_genome_acc}
			 &&
			 $prev_hit->{sorient} eq $curr_hit->{sorient}
			 &&
			 (! &overlap_too_much([$prev_hit_qLEND, $prev_hit_qREND], [$curr_hit_qLEND, $curr_hit_qREND]))
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

	my ($A1, $A2) = @$coordset_A_aref;
	my ($B1, $B2) = @$coordset_B_aref;

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
		return (0);
	}
}

	
####
sub print_chains {
	my @chains = @_;

	foreach my $chain (@chains) {
		print "#\n"; # chain separator
		foreach my $hit (@$chain) {
			print $hit->{line};
		}
	}

	return;
}


		
