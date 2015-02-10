package Fasta::FileChunker;

use strict;
use warnings;
use Carp;

sub new {
	my $packagename = shift;
	my ($outputdirectory) = @_;  ## Full path required.

	unless ($outputdirectory =~ m|^/|) { 
		confess "Error, method requires full path to destination directory";
	}
	
	my $self = {
		outputdirectory => $outputdirectory,
	};

	bless ($self, $packagename);

	unless (-d $outputdirectory) {
		mkdir $outputdirectory or confess "Error, cannot mkdir $outputdirectory";
	}
	
	return ($self);
}


####
sub write_files {
	my $self = shift;
	my ($fasta_reader_obj, $seqsPerChunk, $filenamer_sref) = @_;
	
	unless (ref $filenamer_sref) {
		$filenamer_sref = $self->_get_default_filenamer_routine();
	}
	
	my $index = -1;
	my $fh;
	
	my $output_directory = $self->{outputdirectory};

	my @filenames;

	while (my $seq_obj = $fasta_reader_obj->next()) {
		$index++;
		if ($index % $seqsPerChunk == 0) {
			## start a new file
			close $fh if $fh;
			my $filename = "$output_directory/" . &$filenamer_sref($index+1); # index + 1 because want filenames indexed starting at one instead of zero
			open ($fh, ">$filename") or confess "Error, cannot write to file $filename";
			push (@filenames, $filename);
		}
		my $fasta_format = $seq_obj->get_FASTA_format();
		print $fh $fasta_format;
	}
	close $fh if $fh;

	return (@filenames);
	
}


####
sub clean_output_directory {
	my $self = shift;

	my $output_directory = $self->{outputdirectory};

	my $cmd = "rm -r $output_directory";
	my $ret = system $cmd;
	return ($ret);
}



#########################
### Private Methods #####
#########################

####
sub _get_default_filenamer_routine {
	my $self = shift;

	my $routine = sub { 
		my ($index) = @_;
		
		return ("chunk_$index.fasta");
	};

	return ($routine);
}

1;


	
	
	
