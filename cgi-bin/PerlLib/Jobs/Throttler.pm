package Jobs::Throttler;

use strict;
use warnings;
use Carp;

use Jobs::Runner;

####
sub new {
	my $packagename = shift;
	my ($num_simultaneous_jobs) = @_;  ## area to keep track of individual jobs.  Must exist already.
	## If admin_dir is not set, a default area is chosen under /tmp
	
	my $admin_dir = "/tmp/ThrottlerAdmin$$." . time();
	
	unless ($num_simultaneous_jobs =~ /^\d+$/ && $num_simultaneous_jobs > 0) {
		confess "Error, set num_simultaneous_jobs to a positive integer";
	}
	
	my $self = {
		num_simultaneous_jobs => $num_simultaneous_jobs,
		admin_dir => $admin_dir,
		job_tracker => {},      # maintains list of active jobs by process IDs.
		failed_jobs => [],      # retain jobs that fail.
	};

	bless ($self, $packagename);
  
	return ($self);
}


####
sub launch {
	my $self = shift;
	my @jobs = @_;
	
	## prepare area:
	my $admin_dir = $self->{admin_dir};
	mkdir $admin_dir or confess "Error, cannot mkdir $admin_dir";
	
	my $num_simultaneous_jobs = $self->{num_simultaneous_jobs};
		
	my $num_active_jobs = 0;
	my %active_jobs; # child PID to job being processed.
	
	while (@jobs) {

		while ($num_active_jobs < $num_simultaneous_jobs && @jobs) {
						
			my $next_job = shift @jobs;
			
			my $pid = fork();
			sleep(1); # one-second delay for both child and parent.  Spaces out the children a bit.
			if ($pid) {
				## parent.
				$num_active_jobs++;
				$active_jobs{$pid} = $next_job;
			}
			else {
				## child
				my $job_runner = new Jobs::Runner($next_job, $admin_dir, "ret");
				$next_job->notify_job_start();
				my $ret = $job_runner->run();
				exit($ret);
			}
		}

		my $child_pid = wait();
		if ($child_pid > 0) {
			my $admin_file = "$admin_dir/$child_pid.ret";
			my $job = $active_jobs{$child_pid};
			$self->_process_finished_child($child_pid, $admin_file, $job);
			delete $active_jobs{$child_pid};
			$num_active_jobs--;
		}
	}

	while ( (my $child_pid = wait()) > 0) {
		my $admin_file = "$admin_dir/$child_pid.ret";                                                                                                         
		my $job = $active_jobs{$child_pid} or confess "Error, cannot recover job for pid $child_pid";                                                                                                                   
		$self->_process_finished_child($child_pid, $admin_file, $job);  
	}


	## Cleanup.
	my $ret = system "rm -rf $admin_dir";
	if ($ret) {
		print STDERR "WARNING: could not full cleanup throttler admin directory: $admin_dir\n";
	}
	
	return;
}

####
sub get_failed_jobs {
	my $self = shift;
	
	my @failed_jobs = @{$self->{failed_jobs}};

	return (@failed_jobs);
}


###############################
# Private Methods
###############################


####
sub _process_finished_child {
	my $self = shift;
	my ($child_pid, $admin_file, $job) = @_;

    if (-e $admin_file) {
		open (my $fh, "$admin_file") or confess "Error, cannot open file $admin_file";
		my $exit_val = <$fh>;
		close $fh;
		
		$job->set_exit_status($exit_val);
		if ($exit_val != 0) {
			$self->_add_failed_job($job);
		}
		
		unlink($admin_file); # cleanup.
		
	}
	else {
		$job->set_exit_status(-1);
		$self->_add_failed_job($job);
	}

	$job->notify_job_end();

	return;
}

####
sub _add_failed_job {
	my $self = shift;
	my ($job) = @_;
	
	push (@{$self->{failed_jobs}}, $job);

	return;
}




1; #EOM
	
	
