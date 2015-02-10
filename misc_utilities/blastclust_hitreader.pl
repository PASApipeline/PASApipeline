#!/usr/bin/env perl

use strict;
use warnings;

# originally written by Bernard Suh, tweaked by bhaas

my $record;
my $IDtype;
my $ListSz;
my $SeqID;
my $SeqLen;
my $length;
my $i = 1;
my ($N1, $N2, $HSPL1, $HSPL2, $Score, $PercId);

my $file = $ARGV[0] or die "usage: $0 hits.file.binary\n\n";

my @seq_names;

open(FILE, "<$file");

# 1-byte boolean IDtype -- 1 if numeric IDs; 0 if string IDs
read(FILE, $record, 1);
$IDtype = my_unpack("H2", $record);

# 3-bytes of junk?
read(FILE, $record, 3);

# 4-byte int -- # of seqs if numeric IDs; length of IDs if string IDs
read(FILE, $record, 4);
$ListSz = my_unpack("H2 H2 H2 H2", $record);
if ($IDtype) { # numeric IDs
	for ($i=1; $i <= $ListSz; $i++) {
		read(FILE, $record, 4);
		$SeqID = my_unpack("H2 H2 H2 H2", $record); 
		
	}
	for ($i=1; $i <= $ListSz; $i++) {
		read(FILE, $record, 4);
		$length = my_unpack("H2 H2 H2 H2", $record);   
	}
} else { # string IDs
	read(FILE, $record, $ListSz);
	my @temp = split(" ", unpack("A*", $record));
	@seq_names = @temp;
	
	my $count = scalar @temp;
	for ($i=1; $i <= $count; $i++) {
		read(FILE, $record, 4);
		$length = my_unpack("H2 H2 H2 H2", $record);   
	}
}




until( eof(FILE) ) {
	read(FILE, $record, 4);
	$N1 = my_unpack("H2 H2 H2 H2", $record);
	read(FILE, $record, 4);
	$N2 = my_unpack("H2 H2 H2 H2", $record);
	
	$N1 = $seq_names[$N1];
	$N2 = $seq_names[$N2];
	
	read(FILE, $record, 4);
	$HSPL1 = my_unpack("H2 H2 H2 H2", $record);
	read(FILE, $record, 4);
	$HSPL2 = my_unpack("H2 H2 H2 H2", $record);
	
	read(FILE, $record, 8);
	$Score = unpack("F", $record);
	read(FILE, $record, 8);
	$PercId = unpack("F", $record);
	print "$N1\t$N2\t$HSPL1\t$HSPL2\t$Score\t$PercId\n";
}
close(FILE);

sub my_unpack {
    my $template = shift;
    my $expr = shift;
    return hex join "", reverse unpack($template, $expr);
}

__END__

Output Format

number of 1st sequence (refer to your multi-fasta file where the first seq is 0)
number of 2nd sequence
HSP length of the 1st sequence
HSP length of the 2nd sequence
BLAST score
% identity


Documentation from NCBI:

Format of the hit-list file.

The hit-list file consists of the following parts:

 - header
 - sequence ID list
 - sequence length list
 - hit list

The byte-by-byte layout is platform-dependent; field sizes given here
are true for most UNIX platforms.

A.1. Header.

        4-byte integer  IDtype  1 if numeric IDs; 0 if string IDs
        4-byte integer  ListSz  size of the ID list; if IDs are numeric this
                                is the number of SeqID records, otherwise this
                                is the length of the ID list (in bytes)

A.2. Sequence ID list.

If IDtype is 1 (numeric IDs) then the list is ListSz records of

        4-byte integer  SeqID   sequence ID (numeric)

If IDtype is 0 (string IDs) then the list is a list of records of

        var-length char SeqID   sequence ID (string)
        space (' ')             separator

(total length is ListSz bytes; the number of sequences is equal to the number
of spaces).

A.3. Sequence length list.

This is a list of

        4-byte integer  SeqLen  sequence length

A.4. Hit list.

The list consists of the following records going to the end of file:

        4-byte integer  N1      ordinal number of the 1st sequence
        4-byte integer  N2      ordinal number of the 2nd sequence
        4-byte integer  HSPL1   HSP length on the 1st sequence
        4-byte integer  HSPL2   HSP length on the 2nd sequence
        8-byte float    Score   BLAST score
        8-byte float    PercId  Percent of identical residues
