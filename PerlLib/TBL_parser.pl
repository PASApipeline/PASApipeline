package TBL_parser;

use strict;
use warnings;
use Carp;

####
sub new {
	my $packagename = shift;
	my ($file) = @_;

	my $self = {
		
		file => $file,
		features => [],
	};

	bless ($self, $packagename);

	my @features = $self->_parse_tbl_file($file);
	
	$self->{features} = \@features;
	
	return($self);
}

####
sub get_features {
	my $self = shift;

	return(@{$self->{features}});
}

####
sub _parse_tbl_file {
	my $self = shift;
	my ($file) = @_;

	my $contig;
	my $feature;
	my @features;

	open (my $fh, $file) or confess "Error, cannot open file $file";
	while (<$fh>) {
		#print;
		unless (/\w/) { next; }
		
		chomp;
		if (/^>Feature\s+(.*)$/) {
			$contig = $1;
		}
		else {
			unless (defined $contig) {
				confess "Error, contig name not parsed from tbl header";
			}
			
			my @x = split(/\t/);
			if (defined($x[2]) && $x[2] =~ /\w/) {
				# got a new feature
				my ($end5, $end3) = ($x[0], $x[1]);
				$end5 =~ s/[\<\>]//g;
				$end3 =~ s/[\<\>]//g;
				my $feat_type = $x[2];
				$feature = TBL_feature->new($contig, $feat_type, $end5, $end3);
				push (@features, $feature);
			}
			elsif ($x[0] =~ /\d/ && $x[1] =~ /\d/) {
				$feature->add_end5_end3($x[0], $x[1]);
			}
			elsif ($x[3] && $x[4]) {
				$feature->add_attribute($x[3], $x[4]);
			}
			elsif ($x[3] eq 'pseudo') {
				$feature->{pseudogene} = 1;
			}
			else {
				print STDERR "Warning, cannot parse line $_\n";
			}
		}
	}
	close $fh;
	
	return(@features);
}



package TBL_feature;

use strict;
use warnings;
use Carp;
use Nuc_translator;

sub new {
	my $packagename = shift;
	
	my ($contig_acc, $feat_type, $end5, $end3) = @_;


	
	unless ($contig_acc && $feat_type && $end5 =~ /\d/ && $end3 =~ /\d/) {
		confess "error w/ params";
	}

	
	my $self = {
		contig_acc => $contig_acc,
		feat_type => $feat_type,
		coords => { $end5 => $end3 },
		strand => undef,
		pseudogene => 0,
		attributes => {},
	};

	bless ($self, $packagename);

	$self->_set_strand_by_coords($end5, $end3);

	return($self);
}


####
sub get_feature_sequence {
	my $self = shift;
	my $genome_ref = shift;
	
	my $feature_seq = "";
	
	foreach my $end5 (sort {$a<=>$b} keys %{$self->{coords}}) {
		my $end3 = $self->{coords}->{$end5};

		my ($lend, $rend) = sort {$a<=>$b} ($end5, $end3);
		
		my $seg_seq = substr($$genome_ref, $lend - 1, $rend - $lend + 1);
		
		$feature_seq .= $seg_seq;
	}

	if ($self->{strand} eq '-') {
		$feature_seq = &reverse_complement($feature_seq);
	}

	return($feature_seq);
}



####
sub toString {
	my $self = shift;
	
	my $text = join("\t", $self->{contig_acc}, 
					$self->{feat_type}) . "\n";
	
	$text  .= "\tstrand($self->{strand})\tpseudogene=$self->{pseudogene}\n";
	
	foreach my $end5 (sort keys %{$self->{coords}}) {
		my $end3 = $self->{coords}->{$end5};
		$text .= "\t$end5-$end3\n";
	}
	
	foreach my $key (sort keys %{$self->{attributes}}) {
		my $vals = join(", ", @{$self->{attributes}->{$key}});
		$text .= "\t$key=[$vals]\n";
	}

	return($text);
}

####
sub build_description_line {
	my $self = shift;
	
	my $text = $self->toString();
	
	$text =~ s/\t/ /g;
	
	$text =~ s/\n/\t/g;
	
	
	return($text);
}	


####
sub add_attribute {
	my $self = shift;
	my ($key, $val) = @_;

	push (@{$self->{attributes}->{$key}}, $val);
	
	return;
}

####
sub _set_strand_by_coords {
	my $self = shift;
	my ($end5, $end3) = @_;
	
	if ($end5 < $end3) {
		$self->{strand} = '+';
	}
	elsif ($end5 > $end3) {
		$self->{strand} = '-';
	}
	
	return;
}


####
sub add_end5_end3 {
	my $self = shift;
	my ($end5, $end3) = @_;

	unless (defined $self->{strand}) {
		$self->set_strand_by_coords($end5, $end3);
	}

	$self->{coords}->{$end5} = $end3;
			
	return;
}

####
sub get_feature_ID {
	my $self = shift;
	
	my $type = $self->{feat_type};
	
	if ($type eq 'gene') {
		return($self->get_single_attribute("locus_tag"));
	}
	elsif ($type eq 'CDS') {
		return($self->get_single_attribute("protein_id"));
	}
	else {
		return(undef);
	}
}

####
sub get_single_attribute {
	my $self = shift;
	my ($attribute_key) = @_;

	my $vals_aref = $self->{attributes}->{$attribute_key};

	if (ref $vals_aref) {
		return($vals_aref->[0]);
	}
	else {
		return(undef);
	}
}
	
1; #EOM

