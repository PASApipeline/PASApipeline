#!/usr/bin/env perl

use strict;
use warnings;
use GD;

use lib ($ENV{EUK_MODULES});
use ColorGradient;

my $usage = "usage: $0 m2fmt_file seqLengthsFile orderedGeneAccsFile\n\n";

my $m2fmt_file = $ARGV[0] or die $usage;
my $seqlengths_file = $ARGV[1] or die $usage;
my $orderedGeneAccsFile = $ARGV[2];

## image settings.
my $TEXT_PANEL_WIDTH = 50;
my $CANVAS_WIDTH = 500;
my $SPACE_BETWEEN_GENES = 25;
my $GENE_HEIGHT = 10;

## image vars:
my %COLORS; # r,g,b -> color
my %ALPHA_COLORS; # r,g,b,alpha  alpha value may range from 0 (opaque) to 127 (transparent)
my $MAX_GENE_LENGTH = 0;

my $ORDERED_MATCHES_FLAG = 0; # if 1, draws from a-> b, but not back from b-> a, so one-sided view.


main: {
	


	my %gene_pairs_to_matches;
	my %genes_of_interest;
	
	my @genes;
	if ($orderedGeneAccsFile) {
		@genes = `cat $orderedGeneAccsFile`;
		foreach my $gene (@genes) { $gene =~ s/\s+//g; }
		foreach my $gene (@genes) {
			$genes_of_interest{$gene} = 1;
		}
		%gene_pairs_to_matches = &parse_matches($m2fmt_file, \%genes_of_interest);
	}
	else {
		%gene_pairs_to_matches = &parse_matches($m2fmt_file);
		@genes = &get_gene_accessions(\%gene_pairs_to_matches);
		foreach my $gene (@genes) {
			$genes_of_interest{$gene} = 1;
		}
		
	}
	
	my $num_genes = scalar (@genes);
	
	my %seqlengths = &parse_seqlengths($seqlengths_file, \%genes_of_interest);


	my $image_height = $GENE_HEIGHT * ($num_genes+2) + $SPACE_BETWEEN_GENES * ($num_genes+2);
	my $image_width = $CANVAS_WIDTH + $TEXT_PANEL_WIDTH + 10;;

	my $image = new GD::Image($image_width, $image_height);
	my $bgcolor = &get_color($image, 255,255,255);
	$image->fill($image_width, $image_height, $bgcolor); #Rectangle(1,1, $image_width, $image_height, $bgcolor);
	
	## draw the genes:
	my %gene_acc_to_YbaseCoord = &draw_genes($image, \@genes, \%seqlengths);

	&draw_matches($image, \%gene_acc_to_YbaseCoord, \%gene_pairs_to_matches);
	
	&add_accession_text($image, \@genes, \%gene_acc_to_YbaseCoord);
	
	print $image->png();

	exit(0);

}

####
sub get_gene_accessions {
	my ($gene_pairs_to_matches_href) = @_;

	my %gene_accs;
	foreach my $gene_pair (keys %$gene_pairs_to_matches_href) {
		my ($geneA, $geneB) = split (/$;/, $gene_pair);
		$gene_accs{$geneA} = 1;
		$gene_accs{$geneB} = 1;
	}

	my @genes = keys %gene_accs;
	return (@genes);
}

####
sub parse_matches {
	my ($m2fmt_file, $accs_restricted_href) = @_;

	my %match_pairs;

	open (my $fh, $m2fmt_file) or die "Error, cannot open file $m2fmt_file";
	while (<$fh>) {
		chomp;
		my @x = split (/\t/);
		my ($accA, $accB, $lendA, $rendA, $lendB, $rendB) = ($x[0],$x[1],$x[17],$x[18],$x[20],$x[21]);

		if ($accA eq $accB) { next; } # no self matches
		
		if ($accs_restricted_href) {
			unless ($accs_restricted_href->{$accA} && $accs_restricted_href->{$accB}) { next; }
		}
		
		if ($ORDERED_MATCHES_FLAG) {
			if ($accA gt $accB) {
				## swap info:
				($accA, $lendA, $rendA, $accB, $lendB, $rendB) = ($accB, $lendB, $rendB, $accA, $lendA, $rendA);
			}
		}
		my $match_token = join ("$;", $accA, $accB);
		push (@{$match_pairs{$match_token}},  [$lendA, $rendA, $lendB, $rendB] );
	}
	close $fh;

	return (%match_pairs);
}


####
sub parse_seqlengths {
	my ($file, $accs_want_href) = @_;

	my %lengths;
	
	open (my $fh, $file) or die "Error, cannot open file $file";
	while (<$fh>) {
		chomp;
		my ($seqlength, $acc, @rest) = split (/\s+/);
		unless ($accs_want_href->{$acc}) { next; }
		$lengths{$acc} = $seqlength;
		if ($seqlength > $MAX_GENE_LENGTH) { 
			$MAX_GENE_LENGTH = $seqlength;
		}
	}
	close $fh;
	
	return (%lengths);
}


####
sub seqCoord_to_imageCoord {
	my ($seq_coord) = @_;

	my $ratio_length = $seq_coord / $MAX_GENE_LENGTH;
	my $offset = int($ratio_length * $CANVAS_WIDTH);
	
	my $image_coord = $TEXT_PANEL_WIDTH + $offset;
	return ($image_coord);
}
					 
####
sub draw_genes {
	my ($image, $genes_aref, $seq_lengths_href)  = @_;

	my $y_coord = 0;
	
	my %gene_acc_to_YbaseCoord;

	my $black = &get_color($image, 0,0,0);

	foreach my $gene (@$genes_aref) {
		$y_coord += $GENE_HEIGHT + $SPACE_BETWEEN_GENES;
		$gene_acc_to_YbaseCoord{$gene} = $y_coord;
		
		my $seqlength = $seq_lengths_href->{$gene} or die "Error, no seqlength for $gene";
		my $x_left = &seqCoord_to_imageCoord(1);
		my $x_right = &seqCoord_to_imageCoord($seqlength);
		my $y_coord2 = $y_coord-$GENE_HEIGHT;
		
		$image->filledRectangle($x_left, $y_coord2, $x_right, $y_coord, $black);
		print STDERR "Drawing rect:  ($x_left,$y_coord2) ($x_right,$y_coord)\n";
	}
	
	return (%gene_acc_to_YbaseCoord);
}

####
sub get_color {
	my ($image, $r, $g, $b) = @_;

	my $colortoken = join (",", $r, $g, $b);

	if (my $color = $COLORS{$colortoken}) {
		return ($color);
	}
	else {
		my $color = $image->colorAllocate($r, $g, $b);
		$COLORS{$colortoken} = $color;
		return ($color);
	}
}


####
sub add_accession_text {
	my ($image, $genes_aref, $gene_acc_to_YbaseCoord_href) = @_;
	
	my $color = &get_color($image, 0, 0, 0);

	foreach my $gene (@$genes_aref) {
		my $ybase_coord = $gene_acc_to_YbaseCoord_href->{$gene};
		my $ycoord = $ybase_coord - $GENE_HEIGHT;
		
		$image->string(gdSmallFont, 2, $ycoord, $gene, $color);
	}

	return;
}



####
sub get_alpha_color {
	my ($image, $r, $g, $b, $alpha) = @_;

	my $colortoken = join (",", $r, $g, $b, $alpha);

	if (my $color = $ALPHA_COLORS{$colortoken}) {
		return ($color);
	}
	else {
		my $color = $image->colorAllocateAlpha($r, $g, $b, $alpha);
		$ALPHA_COLORS{$colortoken} = $color;
		return ($color);
	}
}


####
sub draw_matches {
	my ($image, $gene_acc_to_YbaseCoord_href, $gene_pairs_to_matches_href) = @_;
	
	my $num_matches = 0;
	my @gene_accs = (keys %$gene_acc_to_YbaseCoord_href);
	my $num_genes = scalar (@gene_accs);
	my @colors = &ColorGradient::get_RGB_gradient($num_genes);
	
	my %gene_acc_to_color;
	foreach my $gene (@gene_accs) {
		$gene_acc_to_color{$gene} = shift @colors;
	}
	

	foreach my $gene_pair (keys %$gene_pairs_to_matches_href) {
		$num_matches += scalar (@{$gene_pairs_to_matches_href->{$gene_pair}});
	}
	

	foreach my $gene_pair (keys %$gene_pairs_to_matches_href) {
		my ($geneA, $geneB) = split (/$;/, $gene_pair);
		
		my ($r,$g,$b) = @{$gene_acc_to_color{$geneA}};
		my $color = &get_color($image, $r,$g,$b);
		
		my $Ycoord_geneA = $gene_acc_to_YbaseCoord_href->{$geneA};
		my $Ycoord_geneB = $gene_acc_to_YbaseCoord_href->{$geneB};
		if ($Ycoord_geneA > $Ycoord_geneB) {
			$Ycoord_geneA -= $GENE_HEIGHT;
		}
		else {
			$Ycoord_geneB -= $GENE_HEIGHT;
		}
		
		

		
		my @matches = @{$gene_pairs_to_matches_href->{$gene_pair}};
	
		my $num_matches = scalar(@matches);

		

		foreach my $match (@matches) {
			my ($lendA, $rendA, $lendB, $rendB) = @$match;
			
			print STDERR "-processing match coords: @$match\n";
			
			my $xcoord_lendA = &seqCoord_to_imageCoord($lendA);
			my $xcoord_rendA = &seqCoord_to_imageCoord($rendA);
			
			my $xcoord_lendB = &seqCoord_to_imageCoord($lendB);
			my $xcoord_rendB = &seqCoord_to_imageCoord($rendB);
			
			my $pixels_lengthA = $xcoord_rendA - $xcoord_lendA;
			my $pixels_lengthB = $xcoord_rendB - $xcoord_lendB;

			## draw lines:
			my $pixel_step = int(rand(9)) + 3;
			for (my $i = $xcoord_lendA; $i <= $xcoord_rendA; $i+=$pixel_step) {
				
				my $ratio_lenA = ($i - $xcoord_lendA) / $pixels_lengthA;
				my $rel_B = int($ratio_lenA * $pixels_lengthB + 0.5);
				
				my $xcoordB = $xcoord_lendB + $rel_B;

				$image->line($i, $Ycoord_geneA, $xcoordB, $Ycoord_geneB, $color);
			}
			
				

			
			#my $poly = new GD::Polygon;
			#$poly->addPt($xcoord_lendA, $Ycoord_geneA);
			#$poly->addPt($xcoord_rendA, $Ycoord_geneA);
			
			#$poly->addPt($xcoord_lendB, $Ycoord_geneB);
			#$poly->addPt($xcoord_rendB, $Ycoord_geneB);
			
			#$image->filledPolygon($poly, $color);
		}
	}
	
	return;
}
