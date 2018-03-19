package MALIGN::PAML::yn00_output_parser;

use strict;
use warnings;
use Carp;

####
sub get_omega_lines {
	my ($full_yn_text) = @_;
	
	unless ($full_yn_text =~ /\w/) {
		confess "Error, input 'yn' text is lacking";
	}
	

	my @lines = split (/\n/, $full_yn_text);
	while (@lines) {
		my $line = shift @lines;
		if ($line =~ /omega/) {
			shift @lines;
			last;
		}
	}
	
	my @omega_line_structs;
	
	my $line = shift @lines;
	if ($line =~ /\w/) {
		while ($line =~ /\w/) {
			$line =~ s/^\s+|\s+$//g; # trim leading/trailing ws
			my @x = split (/\s+/, $line);
			
			unless (scalar @x == 13) {
				confess "Error, cannot parse omega line: $line ";
			}
			my $omega_struct = { 
				line => $line,
				seq1 => $x[0],
				seq2 => $x[1],
				S => $x[2],
				N => $x[3],
				t => $x[4],
				kappa => $x[5],
				omega => $x[6],
				dN => $x[7],
				dN_SE => $x[9],
				dS => $x[10],
				dS_SE => $x[12]
								 };
			
			push (@omega_line_structs, $omega_struct);
			$line = shift @lines;
		}
		
		return (@omega_line_structs);
	}
	else {
		return;
	}
}

1;
