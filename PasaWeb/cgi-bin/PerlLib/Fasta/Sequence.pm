package Fasta::Sequence;
use strict;

sub new {
    my ($packagename, $header, $sequence) = @_;
    
    ## extract an accession from the header:
    my ($acc, $rest) = split (/\s+/, $header, 2);
        
    my $self = { accession => $acc,
		 header => $header,
		 sequence => $sequence,
		 filename => undef };
    bless ($self, $packagename);
    return ($self);
}

####
sub get_accession {
    my $self = shift;
    return ($self->{accession});
}

####
sub get_header {
    my $self = shift;
    return ($self->{header});
}

####
sub get_sequence {
    my $self = shift;
    return ($self->{sequence});
}

#### 
sub get_FASTA_format {
    my $self = shift;
    my $header = $self->get_header();
    my $sequence = $self->get_sequence();
    $sequence =~ s/(\S{60})/$1\n/g;
    my $fasta_entry = ">$header\n$sequence\n";
    return ($fasta_entry);
}


####
sub write_fasta_file {
    my $self = shift;
    my $filename = shift;

    my ($accession, $header, $sequence) = ($self->{accession}, $self->{header}, $self->{sequence});
    $sequence =~ s/(\S{60})/$1\n/g;
   
    my $tempfile;
    if ($filename) {
	$tempfile = $filename;
    } else {
	my $acc = $accession;
	$acc =~ s/\W/_/g;
	$tempfile = "$acc.fasta";
    }
    
    open (TMP, ">$tempfile") or die "ERROR! Couldn't write a temporary file in current directory.\n";
    print TMP ">$header\n$sequence";
    close TMP;
    return ($tempfile);
}

1; #EOM


