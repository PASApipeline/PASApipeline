package IniReader;


use strict;
use warnings;
use Carp;



sub new {
  my $packagename = shift;
  my ($filename) = @_;

  my $self = { section_to_att_val => {},  # section -> att = value
			   
			   };
  


  open (my $fh, $filename) or confess "Error, cannot open file $filename";
  
  my $current_section = "";

  while (<$fh>) {
	if (/^[\:\#]/) { next; } ## comment line
	unless (/\w/) { next; }
	
	if (/\[([^\]]+)\]/) {
	  $current_section = $1;
	  $current_section = &_trim_flank_ws($current_section);
	}
	elsif (/^(.*)=(.*)$/) {
	  my $att = $1;
	  my $val = $2;

	  $att = &_trim_flank_ws($att);
	  $val = &_trim_flank_ws($val);
	  $self->{section_to_att_val}->{$current_section}->{$att} = $val;
	  
	}
  }
  close $fh;

  bless ($self, $packagename);


  return($self);
}

####
sub get_section_headings {
  my $self = shift;
  my @section_headings = keys %{$self->{section_to_att_val}};

  return(@section_headings);
}


####
sub get_section_attributes {
  my $self = shift;
  my $section = shift;
  
  my @attributes = keys %{$self->{section_to_att_val}->{$section}};
  
  return(@attributes);
}

####
sub get_value {
  my $self = shift;
  my ($section, $attribute) = @_;


  return ($self->{section_to_att_val}->{$section}->{$attribute});
}


####
sub _trim_flank_ws {
  my ($string) = @_;

  $string =~ s/^\s+|\s+$//g;
  
  return($string);
}


1; #EOM

	  

	
