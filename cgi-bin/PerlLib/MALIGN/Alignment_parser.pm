package MALIGN::Alignment_parser;

use strict;
use warnings;
use MALIGN::Alignment;

####
sub parse_alignment {
	my ($alignment_file, $alignment_format) = @_;

	my %supported_formats = ( 'phylip' => 1,
							  'stockholm' => 1,
							  );
	
	unless ($supported_formats{$alignment_format}) {
		die "Error, do not support format $alignment_format, only: " . join (", ", keys %supported_formats);
	}
	
	
	my $alignment_obj = undef;

	if ($alignment_format eq 'phylip') { 
		$alignment_obj = &parse_phylip_alignment($alignment_file);
	}

	elsif ($alignment_format eq 'stockholm') {
		$alignment_obj = &parse_stockholm_alignment($alignment_file);
	}

	
	
	else {
		# shouldn't get here if supported formats check works above.
		die "Error, do not recognize alignment format: $alignment_format\n";
	}

   	return ($alignment_obj);
}


####
sub parse_phylip_alignment {
	my ($alignment_file) = @_;

	open (my $fh, $alignment_file) or die "Error, cannot open file $alignment_file";
	my $header_line = <$fh>;
	$header_line =~ s/^\s+//;
	my ($num_alns, $num_chars) = split (/\s+/, $header_line);
	unless ($num_alns =~ /^\d+$/ && $num_chars =~ /^\d+$/) {
		die "Error, couldn't parse header of phylip formatted file: $header_line";
	}

	my @alignment_char_strings;
	my @accs;
	

	## Read first phylip block
	my $counter = 0;
	while (<$fh>) {
		chomp;
		unless (/\S/) { last; }
		my ($acc, $alignment_string) = split (/\s+/, $_, 2);
		unless ($acc && $alignment_string) {
			die "Error reading phylip alignment block. $_";
		}
		# initialize alignment info
		$accs[$counter] = $acc;
		$alignment_char_strings[$counter] = $alignment_string;
		
		$counter++;
	}
	
	## Read all other phylip blocks:
	$counter = 0;
	while (<$fh>) {
		chomp;
		s/^\s+//; # trim leading whitespace
		if (/\S/) {
			$alignment_char_strings[$counter] .= $_;
			$counter++;
		}
		else {
			$counter = 0; # reinitialize counter for next alignment block;
		}
	}
	close $fh;
	
	## instantiate an alignment object

	my $alignment_obj = MALIGN::Alignment->new();
	
	for (my $i = 0; $i <= $#accs; $i++) {
		my $acc = $accs[$i];
		my $alignment_string = $alignment_char_strings[$i];
		$alignment_string =~ s/\s//g;

		my $aligned_sequence_obj = MALIGN::Aligned_sequence->new($acc, $alignment_string);
		$alignment_obj->add_aligned_sequence($aligned_sequence_obj);
	}


	return ($alignment_obj);
}
	
####
sub parse_stockholm_alignment {
	my  ($alignment_file) = @_;


	my $alignment_obj = MALIGN::Alignment->new();

	open (my $fh, $alignment_file) or die "Error, cannot open file $alignment_file";
	while (<$fh>) {
		unless (/\w/) { next; }
		chomp;
		my ($acc, $aligned_seq) = split (/\s+/);

		$acc =~ s|/\d+-\d+$||; #remove any range indicator, want just the accession.
		$aligned_seq =~ s/\./-/g; #make dash characters are canonical gap representation.


		my $aligned_sequence_obj = new MALIGN::Aligned_sequence($acc, $aligned_seq);
		
		$alignment_obj->add_aligned_sequence($aligned_sequence_obj);

	}
	close $fh;
	
	return ($alignment_obj);
}
		



1;
