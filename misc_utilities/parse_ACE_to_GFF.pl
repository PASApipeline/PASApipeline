#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\nusage: $0 ace_file\n\n";

my $ace_file = $ARGV[0] or die $usage;

main: {
	
	## parse one contig record at a time:
	my $contig_record = "";
		
	open (my $fh, $ace_file) or die "Error, cannot open file $ace_file";
	
	<$fh>; 
	<$fh>; # rid headers before first contig record.

	while (<$fh>) {
		if (/^CO /) {
			# start new contig record:
			&parse_contig_record($contig_record) if $contig_record =~ /^CO /;
			$contig_record = "";
		}
		$contig_record .= $_;
	}
	close $fh;
	&parse_contig_record($contig_record) if $contig_record =~ /^CO /;


	exit(0);
}

####
sub parse_contig_record {
	my ($contig_record) = @_;
	
	## contig record 
	my @sections = split (/\n\n+/, $contig_record);
		
	my $contig_section = shift @sections;

	my ($contig_header, $contig_sequence) = &parse_contig_section($contig_section);
	my @x = split (/\s+/, $contig_header);
	my $contig_acc = $x[1];
	my $num_reads = $x[3];
	my $contig_orient = $x[5];
	unless ($contig_orient eq 'U') { die "Error, unexpected contig orient: $contig_orient, expecting U here.";}
		
	## BQ record
	my $BQ_record = shift @sections;
	# -not doing anything with the BQ record

	## AF and BS sequence position records
	my $seq_position_section = shift @sections;
	
	my %read_info = &parse_read_positions($seq_position_section);
	
	my $num_reads_parsed = scalar (keys %read_info);
	if ($num_reads != $num_reads_parsed) {
		die "Error, should have $num_reads reads, but only parsed $num_reads_parsed "; 
	}
	
	my %padded_read_sequences = &parse_padded_read_sequences(@sections);
	my $num_seqs_parsed = scalar (keys %padded_read_sequences);
	if ($num_seqs_parsed != $num_reads) {
		die "Error, should have $num_reads read sequences but parsed only $num_seqs_parsed";
	}
	
	&convert_to_GFF({ contig_acc => $contig_acc,
					  contig_seq => $contig_sequence,
					  read_positions => \%read_info,
					  read_sequences => \%padded_read_sequences,
				  } );
	
	return;
}


####
sub convert_to_GFF {
	my ($contig_data_href) = @_;
	
	# unwrap
	my $contig_acc = $contig_data_href->{contig_acc};
	my $contig_seq = $contig_data_href->{contig_seq};
	my $contig_length = length($contig_seq);
	my %read_positions = %{$contig_data_href->{read_positions}};
	my %read_seqs = %{$contig_data_href->{read_sequences}};
	

	## compute read end position based on sequence length
	foreach my $read_acc (keys %read_seqs) {
		my $read_seq = $read_seqs{$read_acc};
		my $read_seq_len = length($read_seq);
		$read_positions{$read_acc}->{read_length} = $read_seq_len;
		my $padded_start_pos = $read_positions{$read_acc}->{padded_start};
		$read_positions{$read_acc}->{padded_end} = $padded_start_pos + $read_seq_len - 1;
	}


	my @illustration_structs = ([$contig_acc, 1, $contig_length, '+']);


	## print contig record:
	print join ("\t", $contig_acc, "CAP3", "contig", 1, $contig_length, ".", "+", ".", "$contig_acc") . "\n";
	
	## print the read records:
	foreach my $read_acc (sort {$read_positions{$a}->{padded_start}<=>$read_positions{$b}->{padded_start}} keys %read_positions) {
		my $read = $read_positions{$read_acc};
## unwrap the read:
		my $read_orient = ($read->{orient} eq 'U') ? '+' : '-';
		
		my $start = $read->{padded_start};
		my $end = $read->{padded_end};
		
		print join ("\t", $contig_acc, "CAP3", "read", $start, $end, ".", "$read_orient", ".", "$read_acc") . "\n";
	
		push (@illustration_structs, [$read_acc, $start, $end, $read_orient]);

	}

	print "\n"; ## spacer
	
	&draw_contig_illustration($contig_length, \@illustration_structs);

	return;
}


####
sub draw_contig_illustration {
	my ($contig_length, $illustration_structs_aref) = @_;
	
	my $MAX_ILLUST_CHARS = 60; # canvas width measured as ascii chars.

	foreach my $illustration_struct (@$illustration_structs_aref) {
		my ($acc, $lend, $rend, $orient) = @$illustration_struct;
		my $lend_pos = &compute_ascii_pos($lend, $contig_length, $MAX_ILLUST_CHARS);
		my $rend_pos = &compute_ascii_pos($rend, $contig_length, $MAX_ILLUST_CHARS);
		
		print "#\t" . ' ' x ($lend_pos-1);
		if ($orient eq '+') {
			print "-" x ($rend_pos - $lend_pos -1) , ">" , ' ' x ($MAX_ILLUST_CHARS - $rend_pos) , "\t" , $acc , "\n";
		}
		else {
			print "<" , "-" x ($rend_pos - $lend_pos -1) , ' ' x ($MAX_ILLUST_CHARS - $rend_pos) , "\t" , $acc , "\n"; 
		}
	}
	print "\n"; # spacer
	
	return;
}


####
sub compute_ascii_pos {
	my ($seq_pos, $max_seq_pos, $max_illustration_width) = @_;

	my $ratio = $seq_pos / $max_seq_pos;
	
	my $ascii_pos = int($ratio * $max_illustration_width + 4/9);

	return ($ascii_pos);
}






####
sub parse_padded_read_sequences {
	my @read_sections = @_;

	my %padded_read_seqs;
	foreach my $read_section (@read_sections) {
		if ($read_section =~ /^RD/) {
			my @lines = split (/\n/, $read_section);
			my $read_header = shift @lines;
			my $read_seq = join ("", @lines);
			$read_seq =~ s/\s//g;
			my @x = split (/\s+/, $read_header);
			my $acc = $x[1];
			$padded_read_seqs{$acc} = $read_seq;
		}
	}

	return (%padded_read_seqs);
}


####
sub parse_read_positions {
	my ($seq_position_section) = @_;

	my %read_info;
	
	my @lines = split (/\n/, $seq_position_section);
	foreach my $line (@lines) {
		my @x = split (/\s+/, $line);
		my $token = $x[0];
		if ($token eq 'AF') {
			my ($af, $acc, $orient, $padded_start_pos) = @x;
			$read_info{$acc}->{orient} = $orient;
			$read_info{$acc}->{padded_start} = $padded_start_pos;
		}
		elsif ($token eq 'BS') {
			; # ignoring
		}
		else {
			die "//\n$seq_position_section\n//\nError, don't understand token $token ";
		}
	}

	return (%read_info);
}





####
sub parse_contig_section {
	my ($contig_section) = @_;
		
	my @lines = split (/\n/, $contig_section);
	my $contig_header = shift @lines;
	while ($contig_header !~ /^CO / && @lines) {
		$contig_header = shift @lines;
	}
	unless ($contig_header) { die "//\n$contig_section\n//\nError, didn't parse contig header"; }
	
	## create the contig sequence:
	my $contig_sequence = "";
	while (@lines) {
		my $line = shift @lines;
		unless ($line =~ /\w/) {
			last; # end of contig sequence 
		}
		$contig_sequence .= $line;
	}
	unless ($contig_sequence =~ /\w/) { die "//\n$contig_section\n//\n, Error, didn't recreate contig sequence"; }
	$contig_sequence =~ s/\s//g;
	
	return ($contig_header, $contig_sequence);
}

	
	
	
	
