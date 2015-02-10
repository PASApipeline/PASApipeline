#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Gene_obj;
use Fasta_reader;
use CdbTools;
use Nuc_translator;

my $usage = "usage: $0 genome_seq_lengths_file  six_frame_pep_file genome_seq_file\n\n";

my $genome_seq_lengths_file = $ARGV[0] or die $usage;
my $six_frame_pep_file = $ARGV[1] or die $usage;
my $genome_seq_file = $ARGV[2] or die $usage;

my %genome_lengths;
{ # parse the genome sequence lengths

	open (my $fh, $genome_seq_lengths_file) or die "Error, cannot open $genome_seq_lengths_file";
	while (<$fh>) {
		chomp;
		my ($seq_length, $contig_acc) = split (/\t/);
		$genome_lengths{$contig_acc} = $seq_length;
	}
	close $fh;
}


my %contig_to_gene_lists;


my $fasta_reader = new Fasta_reader($six_frame_pep_file);
while (my $seq_obj = $fasta_reader->next()) {
	
	my $accession = $seq_obj->get_accession();
	my $sequence = $seq_obj->get_sequence();

	print STDERR "\r    parsing $accession           ";

	my ($contig, $pos_info) = split (/\^/, $accession);

	# remove any org-specific prefix from the contig accession
	if ($contig =~ /^(\w\w_)/) {
		my $prefix = $1;
		$contig =~ s/$prefix//;
	}

	my $contig_length = $genome_lengths{$contig} or die "Error, cannot determine length of contig: $contig ";

	my ($strand_frame, $coord) = split (/-/, $pos_info);
	
	my $strand = ($strand_frame =~ /^Fp/) ? '+' : '-';

	$strand_frame =~ /(\d)/;
	my $frame = $1 or die "Error, cannot parse frame from strand_frame: $strand_frame, accession: $accession";
	
	if ($frame < 0 || $frame > 3) {
		die "Error, frame is $frame, unallowed.  acc: $accession";
	}

	my $prot_len = length($sequence);
		
	# base coord on the nucleotide sequence:
	my $start_pos = $coord * 3 + $frame;
	
	my $end_pos = $prot_len * 3 + 3 + $start_pos - 1; # add 3 for stop codon.
	
	if ($end_pos > $contig_length) {
		$end_pos = $contig_length; # stop codon addition might have pushed us over the edge.
	}
	
	if ($strand eq '-') {
		## revcomp the coords:
		$start_pos = $contig_length - $start_pos + 1;
		$end_pos = $contig_length - $end_pos + 1;
	}

	
	#print "$accession\t$contig\t$start_pos\t$end_pos\t$strand\n";

	my ($gene_lend, $gene_rend) = sort {$a<=>$b} ($start_pos, $end_pos);

	push (@{$contig_to_gene_lists{$contig}}, [$accession, $gene_lend, $gene_rend, $strand]);
}

foreach my $contig (sort keys %contig_to_gene_lists) {

	my $genome_seq = &cdbyank_linear($contig, $genome_seq_file);
	
	my @genes = @{$contig_to_gene_lists{$contig}};

	foreach my $gene (@genes) {
		
		my ($accession, $gene_lend, $gene_rend, $orient) = @$gene;
		
		my $gene_seq = substr($genome_seq, $gene_lend - 1, $gene_rend - $gene_lend + 1);

		if ($orient eq '-') {
			$gene_seq = &reverse_complement($gene_seq);
		}

		print ">$accession\n$gene_seq\n";
	}
}




exit(0);
