#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 m2fmt_file  searchDatabase\n\n";

my $m2fmt_file = $ARGV[0] or die $usage;
my $searchDatabase = $ARGV[1] or die $usage;

my %headers;
{
	open (my $fh, $searchDatabase) or die "Error, cannot open file $searchDatabase";
	while (<$fh>) {
		chomp;
		if (/^>(\S+) (.*)/s) {
			my $acc = $1;
			my $description = $2;
			$headers{$acc} = $description;
		}
	}
	close $fh;
}

my $match_counter = 0;
open (my $fh, $m2fmt_file) or die "Error, cannot open file $m2fmt_file";
while (<$fh>) {
	$match_counter++;
	my $match_id = sprintf ("match%05d", $match_counter);
	chomp;
	my @x = split (/\t/);
	my ($contig, $acc, $Evalue, $per_id, $lend, $rend, $orient, $mlend, $mrend) = ($x[0], $x[1], $x[2], $x[10], $x[17], $x[18], $x[16], $x[20], $x[21]);

	($lend, $rend) = sort {$a<=>$b} ($lend, $rend);

	my $strand = ($orient > 0) ? '+' : '-';
	
	my $descr = &get_header($acc);
	$descr =~ s/;/ /g;
	print join ("\t", $contig, "searchDatabase", "translated_nucleotide_match", $lend, $rend, $Evalue, $strand, ".",
				"ID=$match_id; Name=$acc $Evalue $per_id $descr; Target=$acc $mlend $mrend") . "\n";
}
close $fh;


exit(0);
	

####
sub get_header {
	my ($acc) = @_;
	
	if (my $descr = $headers{$acc}) {
		return ($descr);
	}

	else {
		## find acc as substr of existing acc:
		foreach my $stored_acc (keys %headers) {
			if (index($stored_acc, $acc) >= 0) {
				my $descr = $headers{$stored_acc};
				$headers{$acc} = $descr;
				return ($descr);
			}
		}
	}
	
	## if got here, didn't find stored accession.
	die "Error, no description line found for accession: [$acc] ";
}

