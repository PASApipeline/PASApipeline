package Jobs::Job;

use strict;
use warnings;
use Carp;

## Abstract Class, to be inherited from only.


####
sub new {
	my $packagename = shift;

	my $self = { 
		exit_status => undef,
	};
	
	bless ($self, $packagename);

	return ($self);

}

####
sub run {
	my $self = shift;

	confess "Abstract class.\n";

	my $exit_val = 0; ## zero is success, nonzero is failure

	return ($exit_val);
		
}

####
sub set_exit_status {
	my $self = shift;
	my ($exit_status) = @_;

	$self->{exit_status} = $exit_status;
	
	return;
}


####
sub notify_job_start {
	my $self = shift;
	
	## null, placeholder for something good.

	return;
}

####
sub notify_job_end {
	my $self = shift;

	## null, placeholder for something good.

	return;
}



1;

