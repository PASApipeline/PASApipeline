#!/usr/local/bin/perl -w

# lightweight fasta reader capabilities:
package Fasta::Reader;


use strict;
use warnings;
use Fasta::Sequence;

sub new {
    my ($packagename, $fastaFile) = @_;

	## note: fastaFile can be a filename or an IO::Handle
	

    my $self = { fastaFile => undef,,
				 fileHandle => undef };

    bless ($self, $packagename);
    
    ## create filehandle
    my $filehandle = undef;
    
	if (ref $fastaFile eq 'IO::Handle') {
		$filehandle = $fastaFile;
	}
	else {
		
		open ($filehandle, $fastaFile) or die "Error: Couldn't open $fastaFile\n";
		$self->{fastaFile} = $fastaFile;
	}
	
	$self->{fileHandle} = $filehandle;

    return ($self);
}



#### next() fetches next Sequence object.
sub next {
    my $self = shift;
    my $orig_record_sep = $/;
    $/="\n>";
    my $filehandle = $self->{fileHandle};
    my $next_text_input = <$filehandle>;
    my $seqobj = undef;
    if ($next_text_input) {
	$next_text_input =~ s/^>|>$//g; #remove trailing > char.
	$next_text_input =~ tr/\t\n\000-\037\177-\377/\t\n/d; #remove cntrl chars
	my ($header, @seqlines) = split (/\n/, $next_text_input);
	my $sequence = join ("", @seqlines);
	$sequence =~ s/\s//g;
		
	$seqobj = Fasta::Sequence->new($header, $sequence);
    }
    
    $/ = $orig_record_sep; #reset the record separator to original setting.
    
    return ($seqobj); #returns null if not instantiated.
}


#### finish() closes the open filehandle to the query database.
sub finish {
    my $self = shift;
    my $filehandle = $self->{fileHandle};
    close $filehandle;
    $self->{fileHandle} = undef;
}

1; #EOM


