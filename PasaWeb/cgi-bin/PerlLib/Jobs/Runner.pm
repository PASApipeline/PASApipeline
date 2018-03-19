package Jobs::Runner;

use strict;
use warnings;
use Carp;



## Objects of this class execute jobs and write exit values to an administration file in the corresponding admin directory.
##  These admin files are utilized by the Job Throttler, or any other class that wishes to use this utility.

####
sub new {
	my $packagename = shift;
	my ($job_obj, $admin_dir, $file_extension) = @_;

	unless (-d $admin_dir) {
		die "Error, cannot find administration directory: $admin_dir\n";
	}
	
	my $self = {
		job => $job_obj,
		admin_dir => $admin_dir,
		file_ext => $file_extension,
	};
	
	bless ($self, $packagename);

	return ($self);

}

####
sub run {
	my $self = shift;
	
	my $admin_file = $self->_get_admin_file();
	open (my $fh, ">$admin_file") or die "Error, cannot write to $admin_file";
	# run the job.
	my $ret = $self->{job}->run();
	print $fh $ret;
	close $fh;
	
	return ($ret);
}



###################
# Private Methods
###################

sub _get_admin_file {
	my $self = shift;

	my $admin_dir = $self->{admin_dir};
	my $extension = $self->{file_ext};

	my $pid = $$;

	my $admin_file = "$admin_dir/$pid.$extension";

	return ($admin_file);
}


1; #EOM
