#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Gene_obj;
use Fasta_reader;

my $usage = "usage: $0 genome_seq_lengths_file  six_frame_pep_file\n\n";

my $genome_seq_lengths_file = $ARGV[0] or die $usage;
my $six_frame_pep_file = $ARGV[1] or die $usage;


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

my $fasta_reader = new Fasta_reader($six_frame_pep_file);
while (my $seq_obj = $fasta_reader->next()) {
	
	my $accession = $seq_obj->get_accession();
	my $header = $seq_obj->get_header();
	my $sequence = $seq_obj->get_sequence();

	my ($contig, $pos_info) = split (/\^/, $accession);

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

	my $gene_obj = new Gene_obj();
	$gene_obj->populate_gene_object( {$start_pos => $end_pos}, { $start_pos => $end_pos } );
	$gene_obj->{asmbl_id} = $contig;
	$gene_obj->{TU_feat_name} = $accession;
	$gene_obj->{Model_feat_name} = "$accession.model";
	$gene_obj->{com_name} = "$header";
	print $gene_obj->to_GFF3_format() . "\n";
	

}


exit(0);
