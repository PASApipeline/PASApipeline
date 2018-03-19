package MALIGN::Aligned_sequence;
use strict;
use warnings;
use Carp;

sub new {
	my $packagename = shift;
	my ($accession, $sequence) = @_;

	unless ($accession && $sequence) { die "Error, need acc and sequence as constructor params"; }

	my $self = {
		accession => $accession,
		sequence => $sequence,
		length => length($sequence),
	};

	bless ($self, $packagename);
	return ($self);
}

####
sub get_accession {
	my $self = shift;
	return ($self->{accession});
}

####
sub get_sequence {
	my $self = shift;
	return ($self->{sequence});
}

####
sub get_length {
	my $self = shift;
	return ($self->{length});
}


####
sub percent_identity_to {
	my $self = shift;
	my ($other_aligned_seq_obj) = @_;

	if ($self->{length} != $other_aligned_seq_obj->{length}) {
		confess "Trying to compare two aligned sequences of  unequal length.  Are these really from the same alignment?";
	}
	
	my $num_aligned_positions = 0;
	my $num_identical = 0;

	my $self_seq = $self->get_sequence();
	my $other_seq = $other_aligned_seq_obj->get_sequence();

	my @self_chars = split (//, $self_seq);
	my @other_chars = split (//, $other_seq);

	for (my $i = 0; $i < $self->{length}; $i++) {
		
		if ( $self_chars[$i] !~ /[\.\-X]/i || $other_chars[$i] !~ /[\.\-X]/i) {
			# got a relevant char in either sequence
			$num_aligned_positions++;
			
			if ( $self_chars[$i] !~ /[\.\-X]/i  && $self_chars[$i] eq $other_chars[$i]) {
				$num_identical++;
			}

		}
	}
	
	my $percent_identity = $num_identical / $num_aligned_positions * 100;
	
	return ($percent_identity);

}

####
sub convert_to_codon_align {
	my $self = shift;
	my ($cds_sequence) = @_;

	unless ($cds_sequence) { 
		confess "Error, need cds sequence as parameter";
	}
	
	## parse Codon mapping:
		
	my @chars = split (//, $cds_sequence);
                
	my @codons;
	for (my $i = 0; $i <= $#chars-2; $i+=3) {
		my $codon = join ("", @chars[$i..$i+2]);
		push (@codons, $codon);
	}
	
	## convert alignment now
	my $codon_aln_txt = "";
	
	my $aligned_prot_sequence = $self->get_sequence();
	my @alignment_chars = split (//, $aligned_prot_sequence);
	
	my $codon_counter = -1;
	foreach my $alignment_char (@alignment_chars) {
		if ($alignment_char eq '-' || $alignment_char eq '.') {
			$codon_aln_txt .= "-" x 3;
		}
		else {
			$codon_counter++;
			my $codon = $codons[$codon_counter] or die "Error, no codon found at position $codon_counter of @alignment_chars";
			$codon_aln_txt .= $codon;
		}
	}
	
	my $codon_aligned_seq = new MALIGN::Aligned_sequence($self->get_accession(), $codon_aln_txt);

	return ($codon_aligned_seq);
}


	
1;
