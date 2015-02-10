#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");

use Fasta_reader;
use TBL_parser;

my $usage = "usage: $0 genome.tbl genome.fasta\n\n";

my $genome_tbl_file = $ARGV[0] or die $usage;
my $genome_fasta_file = $ARGV[1] or die $usage;


main: {
	
	my $tbl_reader = TBL_parser->new($genome_tbl_file);

	my @features = $tbl_reader->get_features();

	my %genome = &parse_fasta($genome_fasta_file);
	
	my $counter = 0;
	
	foreach my $feature (@features) {
		#print $feature->toString() . "\n";
		
		my $contig = $feature->{contig_acc};
		
		my $genome_seq = $genome{$contig} or confess "Error, no genome sequence for contig: [$contig] ";
		
		my $sequence = $feature->get_feature_sequence(\$genome_seq);
		
		my $feature_ID = $feature->get_feature_ID();
		
		unless ($feature_ID) {
			$feature_ID = ++$counter;
		}

		&print_GFF($feature_ID, $feature);

	}
	
	exit(0);

}
	

sub print_GFF {
	my ($feature_ID, $feature) = @_;

	my $feat_type = $feature->{feat_type};
	
	# unless ($feat_type eq 'CDS') { return; }

	my $contig_acc = $feature->{contig_acc};
	
	my $strand = $feature->{strand};
	
	my %coords = %{$feature->{coords}};

	my $descr_txt = "";
	foreach my $key (sort keys %{$feature->{attributes}}) {
		my $vals = join(", ", @{$feature->{attributes}->{$key}});
		$descr_txt .= "\t$key=[$vals]\n";
	}

	$descr_txt =~ s/\s+/ /g;
	
	foreach my $end5 (sort {$a<=>$b} keys %coords) {
		my $end3 = $coords{$end5};
		
		my ($lend, $rend) = sort {$a<=>$b} ($end5, $end3);
		
		my $feat_val = ($feat_type eq 'CDS') ? 'mRNA' : $feat_type;
		
		print join("\t", $contig_acc, "GenBank", $feat_type, $lend, $rend, ".", $strand, ".", "$feat_val $feature_ID; $descr_txt") . "\n";
	}
	
	return;
}


####
sub parse_fasta {
	my ($fasta_filename) = @_;

	my $fasta_reader = new Fasta_reader($fasta_filename);
	
	my %genome;

	while (my $seq_obj = $fasta_reader->next() ) {

		my $acc = $seq_obj->get_accession();

		my $sequence = $seq_obj->get_sequence();
		
		#my @parts = split(/\|/, $acc);
		#if (scalar(@parts) > 2) {
		#	shift @parts;
		#	shift @parts;
		#	$acc = join("|", @parts) . "|";
		#}
		
		$genome{$acc} = $sequence;
		print STDERR "-storing [$acc] = ". length($sequence) . "\n";
	}

	return(%genome);
}
		



