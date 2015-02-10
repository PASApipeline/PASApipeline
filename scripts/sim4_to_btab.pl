#!/usr/bin/env perl

use strict;

## searching a cDNA (first param) against a genomic database (second param)
## query = cdna, db = genomic.

## MUST use A=4 option!!!!! (first seq coordinates always referrential)

$/ = "\nseq1";
my $data = <STDIN>; #rid first line.
my $chain_num = 0; #alignment chain.

while ($data = <STDIN>) {
    #print $data;
    my $complement = 0;
    if ($data =~ /\(complement\)/) {
	$complement = 1;
    } 
    my @lines = split (/\n/, $data);
    my ($query, $query_length, $query_name) = &get_name_and_length (shift @lines);
    my ($db, $db_length, $db_accession) = &get_name_and_length (shift @lines);
    my $segment_num = 0;
    my $chain_incrementer_flag = 0;
    ## process alignment chain data, if exists.
    while (@lines) {
	my $line = shift (@lines);
	if ($line =~ /(\d+)-(\d+)\s+\((\d+)-(\d+)\)\s+(\d+)\%/) {
	    my ($query1, $query2, $db_end5, $db_end3, $per_id) = ($1, $2, $3, $4, $5);
	    my $match_length = abs($query1 - $query2) + 1;
	    my ($query_end5, $query_end3) = ($query1, $query2);
	    if ($complement) {
		#database sequence always forward orientation and forward coordinates.
		($db_end5, $db_end3) = &revcomp_coord($db_length, $db_end5, $db_end3);

	    }
	    $segment_num++;
	    if (!$chain_incrementer_flag) {
		$chain_num++;
		$chain_incrementer_flag = 1;
	    }
	    
	    ## create btab line.
	    my @btab;
	    $btab[0] = $db_accession;
	    $btab[2] = $query_length;
	    $btab[3] = "sim4";
	    $btab[4] = $db;
	    $btab[5] = ($query_name) ? $query_name : $query;
	    $btab[6] = $db_end5;
	    $btab[7] = $db_end3;
	    $btab[8] = $query_end5;
	    $btab[9] = $query_end3;
	    $btab[10] = $per_id;
	    $btab[13] = $chain_num;
	    $btab[14] = $segment_num;
	    $btab[18] = $match_length;
	    my $btab_line = join ("\t", @btab);
	    print "$btab_line\n";
	}
    }
}

####
sub revcomp_coord {
    my $seq_length = shift;
    my @coords = @_;
    my @ret_coords;
    foreach my $coord (@coords) {
	my $converted_coord = $seq_length - $coord + 1;
	push (@ret_coords, $converted_coord);
    }
    return (@ret_coords);
}


####
sub get_name_and_length {
    my ($line) = @_;
    if ($line =~ /=\s([^,]+),\s(\d+)\sbp/ ) {
	
    my $filename = $1;
    my $length = $2;
    my $accession;
    if ($filename =~ /^(.+)\((.+)\)/) {
	$filename = $1;
	$accession = $2;
    }
    $filename =~ s/.*\///g;
    $accession =~ tr/\(\)//;
    return ($filename, $length, $accession);
} else {
    return ("","");
}
}
