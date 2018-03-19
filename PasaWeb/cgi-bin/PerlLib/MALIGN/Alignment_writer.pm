package MALIGN::Alignment_writer;

use strict;
use warnings;
use Carp;

my %supported_formats = ("phylip" => 1, 
						 "nexus" => 1);

sub get_alignment_text {
	my ($alignment_obj, $format) = @_;
	
	unless ($supported_formats{$format}) {
		confess "Error, format $format is not supported. Only " . join (", ", keys %supported_formats);
	}

	my $alignment_text = "";

	if ($format eq 'phylip') {
		$alignment_text = &get_phylip_text($alignment_obj);
	}

	elsif ($format eq 'nexus') {
		$alignment_text = &get_nexus_text($alignment_obj);
	}
	
	else {
		confess "format not found. ";
		# shouldn't ever get here.
	}

	return ($alignment_text);
}



####
sub get_phylip_text {
	my ($alignment_obj) = @_;
	
	my $phylip_txt = "";
	
	$phylip_txt .= $alignment_obj->get_num_aligned_seqs() . " " .  $alignment_obj->get_num_chars() . "\n";
	
	my @aligned_seqs = $alignment_obj->get_aligned_seqs();
	
	foreach my $aligned_seq (@aligned_seqs) {
		my $accession = $aligned_seq->get_accession();
		my $sequence = $aligned_seq->get_sequence();
		
		$phylip_txt .= sprintf ("%-10s", $accession) . "  $sequence" . "\n";
		
	}
	
	return ($phylip_txt);
}

####
sub get_nexus_text {
	my ($alignment_obj) = @_;
	
	my $nexus_txt = "";
	
	$nexus_txt .= "#NEXUS\n"
		. "Begin data;\n"
		. "     Dimensions ntax=" . $alignment_obj->get_num_entries() . " nchar=" . $alignment_obj->get_num_chars()  .
		";\n"
		. "     Format datatype=DNA gap=-;\n"
		. "   Matrix\n";
	
	my @entries = $alignment_obj->get_entries();
	foreach my $aligned_seq (@entries) {
		my $accession = $aligned_seq->get_accession();                
		my $sequence =  $aligned_seq->get_sequence();
		
		$nexus_txt .= "$accession  $sequence\n";		
	}
	
	return ($nexus_txt);
}


1;
