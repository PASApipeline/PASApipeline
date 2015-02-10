#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Data::Dumper;

my $usage = "usage: $0 alignments.gff3 [min_per_id]\n\n";


my $trans_gff3_file = $ARGV[0] or die $usage;
my $min_per_id = $ARGV[1] || 0;

main: {
	
	my %genome_trans_to_coords;

    my %trans_id_info;
	
	open (my $fh, $trans_gff3_file) or die "Error, cannot open file $trans_gff3_file";
	while (<$fh>) {
		chomp;
		
		unless (/\w/) { next; }
		
		my @x = split(/\t/);

		unless (scalar (@x) >= 8 && $x[8] =~ /ID=/ && $x[8] =~ /(Name|Target)=/) {
			print STDERR "ignoring line: $_\n";
			next;
		}
		
        if ($x[2] eq "expressed_sequence_match") { next; } # ignoring parent feature

		my $scaff = $x[0];
		my $type = $x[2];
		my $lend = $x[3];
		my $rend = $x[4];
        my $per_id = $x[5];
        
		my $orient = $x[6];
		
		my $info = $x[8];
		
		my @parts = split(/;/, $info);
		my %atts;
		foreach my $part (@parts) {
			$part =~ s/^\s+|\s+$//;
			$part =~ s/\"//g;
			my ($att, $val) = split(/=/, $part);
			
			if (exists $atts{$att}) {
				die "Error, already defined attribute $att in $_";
			}
			
			$atts{$att} = $val;
		}

		my $gene_id = $atts{ID} or die "Error, no gene_id at $_, info:[$info] " . Dumper(\%atts);
		if ($atts{Parent}) {
            $gene_id = $atts{Parent};
        }
        
        my $trans_id = $atts{Target} || $atts{Name}  or die "Error, no trans_id at $_, info:[$info] " . Dumper(\%atts);
		{
			my @pieces = split(/\s+/, $trans_id);
			$trans_id = shift @pieces;
		}
		
		my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);

		$genome_trans_to_coords{$scaff}->{$gene_id}->{$trans_id}->{$end5} = $end3;

        my $seg_len = abs($rend - $lend) + 1;
        
        $trans_id_info{$trans_id}->{sum_len} += $seg_len;
        $trans_id_info{$trans_id}->{sum_per_id_len} += ($seg_len * $per_id);
        

	}


	## Output genes in gff3 format:

	foreach my $scaff (sort keys %genome_trans_to_coords) {

		my $genes_href = $genome_trans_to_coords{$scaff};

		foreach my $gene_id (keys %$genes_href) {

			my $trans_href = $genes_href->{$gene_id};

			foreach my $trans_id (keys %$trans_href) {

                my $trans_info_struct = $trans_id_info{$trans_id} or die "Error, no info stored for $trans_id";
                                
                my $avg_per_id = $trans_info_struct->{sum_per_id_len} / $trans_info_struct->{sum_len};
                if ($avg_per_id < $min_per_id) { next; }
                
				my $coords_href = $trans_href->{$trans_id};

				my $gene_obj = new Gene_obj();

				$gene_obj->{TU_feat_name} = $gene_id;
				$gene_obj->{Model_feat_name} = $trans_id;
				$gene_obj->{com_name} = "$gene_id $trans_id";
				
				$gene_obj->{asmbl_id} = $scaff;
				
				$gene_obj->populate_gene_object($coords_href, $coords_href);
			
				print $gene_obj->to_BED_format(score => int($avg_per_id+0.5));
								
			}
		}
	}


	exit(0);
}

