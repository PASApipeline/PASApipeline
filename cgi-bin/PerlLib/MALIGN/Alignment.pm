package MALIGN::Alignment;

use strict;
use warnings;

use Data::Dumper;

sub new {
	my $packagename = shift;
		
	my $self = {
		aligned_seqs => [],
        num_aligned_seqs => 0,
        num_chars => 0,
     };

     bless ($self, $packagename);
     return ($self);
}
 
####
sub get_aligned_seqs {
	my $self = shift;
	return (@{$self->{aligned_seqs}});
}

####
sub get_num_aligned_seqs {
	my $self = shift;
	return ($self->{num_aligned_seqs});
}

####
sub get_num_chars {
	my $self = shift;
	return ($self->{num_chars});
}





####
sub add_aligned_sequence {
	my $self = shift;
	my ($aligned_seq_obj) = @_;
	
	push (@{$self->{aligned_seqs}}, $aligned_seq_obj);
	$self->{num_aligned_seqs}++;
	
	if ($self->{num_chars} == 0) {
		$self->{num_chars} = $aligned_seq_obj->{length};
	}
	elsif ($self->{num_chars} != $aligned_seq_obj->{length}) {
		die Dumper ($self) . "\nError, trying to add an aligned sequence with inconsistent length to others in the alignment";
	}
	

	return;
}



1; #EOM


