#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use BinaryFeatureSearch;
use GFF3_utils;
use Gene_obj;

my $usage = "usage: $0 annotation.gff3 pasa_assemblies.gff3\n\n";

my $annot_gff3 = $ARGV[0] or die $usage;
my $pasa_assemblies_gff3 = $ARGV[1] or die $usage;

main: {

	## these are formatted for binary feature searches
	my %scaffold_strand_to_annot_exons;
	my %scaffold_strand_to_annot_introns;

	&parse_introns_n_exons_from_annot_gff3($annot_gff3, \%scaffold_strand_to_annot_exons, \%scaffold_strand_to_annot_introns);
	
	# simple ordered list of intron coordinates for each stranded-scaffold.
	
	my %scaffold_strand_to_pasa_exons;
	my %scaffold_strand_to_pasa_introns;
	
	&parse_introns_n_exons_from_alignment_gff3($pasa_assemblies_gff3, \%scaffold_strand_to_pasa_exons, \%scaffold_strand_to_pasa_introns);
	
	
	## examine introns found.
	## create simple lookups
	
	my %scaffold_strand_to_annot_intron_map;
	
	
	my %all_introns;
	
	foreach my $scaffold_strand (keys %scaffold_strand_to_annot_introns) {
		my @introns = @{$scaffold_strand_to_annot_introns{$scaffold_strand}};
		foreach my $intron (@introns) {
			my $lend = $intron->{lend};
			my $rend = $intron->{rend};
			my $key = join("$;", $scaffold_strand, $lend, $rend);
			$scaffold_strand_to_annot_intron_map{$key}++;
			$all_introns{$key} = "annot";
		}
	}

	my %scaffold_strand_to_pasa_intron_map;
	foreach my $scaffold_strand (keys %scaffold_strand_to_pasa_introns) {
		my @introns = @{$scaffold_strand_to_pasa_introns{$scaffold_strand}};
		foreach my $intron (@introns) {
			my ($lend, $rend) = ($intron->{lend}, $intron->{rend});
			my $key = join("$;", $scaffold_strand, $lend, $rend);
			$scaffold_strand_to_pasa_intron_map{$key}++;
			if ($all_introns{$key}) {
				$all_introns{$key} = "both";
			}
			else {
				$all_introns{$key} = "pasa";
			}
		}
	}

	my %type_counter;
	foreach my $intron (keys %all_introns) {
		my $type = $all_introns{$intron};
		$type_counter{$type}++;
	}

	## summarize findings thus far:
	my $num_annot_introns = scalar (keys %scaffold_strand_to_annot_intron_map);
	my $num_pasa_introns = scalar(keys %scaffold_strand_to_pasa_intron_map);
	
	my $num_both = $type_counter{"both"};
	
	print "Percent of $num_annot_introns annotated introns found: " . sprintf("%.2f", $num_both/$num_annot_introns*100) . "\n";
	print "Percent of $num_pasa_introns pasa introns annotated: " . sprintf("%.2f", $num_both/$num_pasa_introns*100) . "\n";
	
	
	## Examine intron readthru: annotated introns that correspond to reconstructed PASA transcript exons.
	my $num_introns_readthru = 0;
	my $num_introns_encroached = 0;
	
	foreach my $scaffold_strand (keys %scaffold_strand_to_annot_introns) {
		
		my @annot_introns = @{$scaffold_strand_to_annot_introns{$scaffold_strand}};
		
		unless (exists $scaffold_strand_to_pasa_exons{$scaffold_strand}) { next; } 

		my @pasa_exons = @{$scaffold_strand_to_pasa_exons{$scaffold_strand}};

		my $readthru_flag = 0;
		my $encroachment_flag = 0;

		foreach my $annot_intron (@annot_introns) {
			
			my ($intron_lend, $intron_rend) = ($annot_intron->{lend}, $annot_intron->{rend});
			
			foreach my $pasa_exon (@pasa_exons) {
				my ($exon_lend, $exon_rend) = ($pasa_exon->{lend}, $pasa_exon->{rend});

				if ($exon_lend > $intron_rend) { last; }
				
				
				## check for overlap
				if ($intron_lend <= $exon_rend && $intron_rend >= $exon_lend) {
					## overlap
					$encroachment_flag++;
					
					if ($exon_lend <= $intron_lend && $intron_rend <= $exon_rend) {
						## containment
						$readthru_flag++;
						last;
					}
				}
			}
		}
		if ($readthru_flag) {
			$num_introns_readthru++;
		}
		elsif ($encroachment_flag) {
			$num_introns_encroached++;
		}
		
	}

	print "$num_introns_readthru of $num_annot_introns (" . sprintf("%.2f", $num_introns_readthru/$num_annot_introns*100) . "%) of annotated introns are readthru by pasa assemblies.\n";

	print "$num_introns_encroached of $num_annot_introns (" . sprintf("%.2f", $num_introns_encroached/$num_annot_introns*100) . "%) of annotated introns are encroached by pasa assemblies.\n";
	


	exit(0);
}


####
sub parse_introns_n_exons_from_annot_gff3 {
	my ($annot_gff3, $exons_href, $introns_href) = @_;

	my $gene_obj_indexer_href = {};
	
	my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($annot_gff3, $gene_obj_indexer_href);

	foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
		
		my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};

		foreach my $gene_id (@gene_ids) {

			my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};

			my $orient = $gene_obj_ref->get_orientation();

			my $stranded_scaffold = "$asmbl_id$orient";

			my @introns = $gene_obj_ref->get_intron_coordinates();
			
			my $gene_id = $gene_obj_ref->{TU_feat_name};

			foreach my $intron (@introns) {
				my ($lend, $rend) = sort {$a<=>$b} @$intron;
				
				push (@{$introns_href->{$stranded_scaffold}}, { acc => $gene_id,
															  lend => $lend,
															  rend => $rend, 
					  } );
			}

			my @exons = $gene_obj_ref->get_exons();
			foreach my $exon (@exons) {
				my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();

				push (@{$exons_href->{$stranded_scaffold}}, { acc => $gene_id,
															lend => $lend,
															rend => $rend, 
					  } );
			}
		}
	}

	## sort everything
	foreach my $list_aref (values %$introns_href) {
		@$list_aref = sort {$a->{lend}<=>$b->{lend}} @$list_aref;
	}
	
	foreach my $list_aref (values %$exons_href) {
		@$list_aref = sort {$a->{lend}<=>$b->{lend}} @$list_aref;
	}
	
	

	return;
}


####
sub parse_introns_n_exons_from_alignment_gff3 {
	my ($pasa_assemblies_gff3, $pasa_exons_href, $pasa_introns_href) = @_;
	
	my %scaffold_strand_transcript_to_exons;
	open (my $fh, $pasa_assemblies_gff3) or die "Error, cannot open file $pasa_assemblies_gff3";
	while (<$fh>) {
		chomp;
		my @x = split(/\t/);
		my $acc_info = $x[8];
		
		$acc_info =~ /ID=([^;]+)/ or die "Error, cannot parse ID from $acc_info";
		my $acc = $1 or die "Error, no acc from $acc_info";

		my $orient = $x[6];
		my $lend = $x[3];
		my $rend = $x[4];

		my $scaffold = $x[0];
		
		my $scaffold_strand_transcript = join("$;", $scaffold, $orient, $acc);
		
		push (@{$scaffold_strand_transcript_to_exons{$scaffold_strand_transcript}}, [$lend, $rend]);
		push (@{$pasa_exons_href->{"$scaffold$orient"}}, { acc => $acc,
																	lend => $lend, 
																	rend => $rend,
			  } );
		
	}
	close $fh;

	my %scaffold_strand_to_introns;
	foreach my $scaff_strand_transcript (keys %scaffold_strand_transcript_to_exons) {
		my ($scaff, $strand, $transcript) = split(/$;/, $scaff_strand_transcript);
		
		my @exons = sort {$a->[0]<=>$b->[0]} @{$scaffold_strand_transcript_to_exons{$scaff_strand_transcript}};

		my $num_exons = scalar(@exons);
		
		if (scalar @exons > 1) {
			
			#print "got $num_exons exons from $scaff $strand $transcript\n";
			
			my $exon = shift @exons;
			while (@exons) {
				my $next_exon = shift @exons;
				my $prev_rend = $exon->[1];
				my $next_lend = $next_exon->[0];

				my $intron_lend = $prev_rend + 1;
				my $intron_rend = $next_lend - 1;

				$scaffold_strand_to_introns{"$scaff$strand"}->{"$intron_lend-$intron_rend"}++;
						
				$exon = $next_exon;
			}
		}
	}

	## Reformat data structure to return.
	
	my $counter = 0;
	foreach my $scaff_strand (keys %scaffold_strand_to_introns) {
		my @introns = keys %{$scaffold_strand_to_introns{$scaff_strand}};
		
		my @intron_coords;
		foreach my $intron (@introns) {
			$counter++;
			my ($lend, $rend) = split(/-/, $intron);
			push (@intron_coords, {acc => $counter,
								   lend => $lend, 
								   rend => $rend,
				  });
		}
		
		@intron_coords = sort {$a->{lend}<=>$b->{lend}} @intron_coords;
		
		$pasa_introns_href->{$scaff_strand} = \@intron_coords;
	}


	foreach my $list_aref (values %$pasa_exons_href) {
		@$list_aref = sort {$a->{lend}<=>$b->{lend}} @$list_aref;
	}
	
	
	return;
}
							 
