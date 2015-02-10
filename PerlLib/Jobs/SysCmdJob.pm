package Jobs::SysCmdJob;

use strict;
use warnings;
use Carp;

use Jobs::Job;
use base ("Jobs::Job");

####
sub new {
	my $packagename = shift;
	my ($cmd_string) = @_;

	
	
	my $self = $packagename->SUPER::new(@_);

	$self->{cmd_string} = $cmd_string;
	
	$self->{output} = undef;

	return ($self);
}

####
sub run {
	my $self = shift;
	
	my $cmd = $self->{cmd_string};

    my $output = `$cmd`;
	my $ret = $?;

	print $output;
	
	$self->{output} = $output;
	
	return ($ret);
}


####
sub notify_job_start {
	my $self = shift;

	print "Starting cmd: " . $self->{cmd_string} . "\n";
	
	return;
}

####
sub notify_job_end {
	my $self = shift;
	
	print "Job ended: " . $self->{cmd_string} . "\n";
	
	return;
}


1; #EOM

