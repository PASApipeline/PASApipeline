#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\nusage: $0 m2fmt_result_files_list.file\n\n";

my $m2fmt_result_files = $ARGV[0] or die $usage;

my %query_to_best_hit;

main: {
	my @result_files = `cat $m2fmt_result_files`;
	chomp @result_files;

	foreach my $result_file (@result_files) {
		&parse_result_file($result_file);
	}

	foreach my $acc (keys %query_to_best_hit) {
		my $best_org_matches_href = $query_to_best_hit{$acc};
		
		$acc =~ /^([^_]+)/ or die "Error, cannot extract org from acc $acc";
		my $org = $1;

		foreach my $other_org (keys %$best_org_matches_href) {

			my $match_struct = $best_org_matches_href->{$other_org};
			my $other_acc = $match_struct->{acc};
			
			## Check to see what this organisms best hit is:
			my $other_org_best_hit = $query_to_best_hit{$other_acc}->{$org};
			if (ref $other_org_best_hit && $other_org_best_hit->{acc} eq $acc) {
				&report_orthologs($match_struct, $other_org_best_hit);
			}
		}
	}
		

	exit(0);
}
	
####
sub parse_result_file {
	my ($result_file) = @_;
	
	print STDERR "\r-parsing file $result_file                    ";
	open (my $fh, $result_file) or die "Error, cannot open file $result_file";
	while (<$fh>) {
		my @x = split (/\t/);
		my $query_acc = $x[0];
		my $db_acc = $x[1];
		my $E_value = $x[2];
		my $score = $x[4];
		my $per_ID = $x[10];

		if ($query_acc eq $db_acc) { next; }
		
		$query_acc =~ /^([^_]+)/ or die "Error, cannot parse org token from $query_acc";
		my $query_org = $1;
		
		$db_acc =~ /^([^_]+)/ or die "Error, cannot parse org token from $db_acc";
		my $db_org = $1;
		
		if ($query_org eq $db_org) { next; }

		if ( (! exists $query_to_best_hit{$query_acc}->{$db_org}) || ($query_to_best_hit{$query_acc}->{$db_org}->{score} < $score) ) {
			$query_to_best_hit{$query_acc}->{$db_org} = { acc => $db_acc,
														  org => $db_org,
														  Evalue => $E_value,
														  score => $score,
														  per_ID => $per_ID, 
													  };
			
		}
		
		
	}
	close $fh;
	
	return;
}

my %pair_reported;
####
sub report_orthologs {
	my ($structA, $structB) = @_;

	($structA, $structB) = sort {$a->{org} cmp $b->{org}} ($structA, $structB);


	my ($accA, $orgA, $EvalueA, $scoreA, $per_ID_A) = ($structA->{acc},
													   $structA->{org},
													   $structA->{Evalue},
													   $structA->{score},
													   $structA->{per_ID});

	my ($accB, $orgB, $EvalueB, $scoreB, $per_ID_B) = ($structB->{acc},
													   $structB->{org},
													   $structB->{Evalue},
													   $structB->{score},
													   $structB->{per_ID});

	my $pair_token = join (";;;", $accA, $accB);
	if ($pair_reported{$pair_token}) {
		return; # already reported them.
	}
	else {
		$pair_reported{$pair_token}++;
		print "$orgA\t$accA\t$scoreA\t$EvalueA\t$per_ID_A"
			. "\t$orgB\t$accB\t$scoreB\t$EvalueB\t$per_ID_B\n";
	}

	return;
}


													  
