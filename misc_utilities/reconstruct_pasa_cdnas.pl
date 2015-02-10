#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Gene_obj;
use CdbTools;

my $usage = "usage: $0 pasa_gff3 genome_sequence\n\n";

my $pasa_gff = $ARGV[0] or die $usage;
my $genome_seq_fasta = $ARGV[1] or die $usage;

my %data;

open (my $fh, $pasa_gff) or die "Error, cannot open file $pasa_gff";
while (<$fh>) {
  chomp;
  my @x = split (/\t/);
  my ($asmbl_id, $lend, $rend, $orient, $gene_info) = ($x[0], $x[3], $x[4], $x[6], $x[8]);
  
  $gene_info =~ /Target=(\S+)/ or die "Error, cannot parse target from gene info: $gene_info";
  my $acc = $1;

  my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);

  $data{$asmbl_id}->{$acc}->{$end5} = $end3;

}

foreach my $asmbl_id (keys %data) {
  my $genome_seq = &cdbyank_linear($asmbl_id, $genome_seq_fasta);
  
  foreach my $acc (keys %{$data{$asmbl_id}}) {

	my $coords_href = $data{$asmbl_id}->{$acc};
	
	my $gene_obj = new Gene_obj();
	
	$gene_obj->populate_gene_object($coords_href, $coords_href, \$genome_seq);
	my $spliced_seq = $gene_obj->get_CDS_sequence();
	print ">$acc\n$spliced_seq\n";
  	
  }
}


exit(0);



