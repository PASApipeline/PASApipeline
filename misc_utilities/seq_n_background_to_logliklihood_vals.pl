#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fasta_reader;
use Nuc_translator;

## hexamer stats
my %framed_hexamers;
my %background_hexamers;


## pentamer stats
my %framed_pentamers;
my %background_pentamers;
my %framed_all_pentamer_counts;

my $usage = "usage: $0 targetCDSs backgroundSeqs\n\n";

my $target_CDS = $ARGV[0] or die $usage;
my $background_seqs = $ARGV[1] or die $usage;


main: {

	&parse_targetCDSs($target_CDS);

	&parse_background($background_seqs);

	&add_pseudocounts();

	&report_logliklihood_ratios();
	
	exit(0);
}





####
sub report_logliklihood_ratios {
	

	## Markov-based probabilities (5th order markov chain):
	
	foreach my $framed_hexamer (sort keys %framed_hexamers) {
		my ($hexamer, $frame) = split (/-/, $framed_hexamer);
		
		my $pentamer = substr($hexamer, 0, 5);

		my $framed_hexamer_count = $framed_hexamers{$framed_hexamer};
		my $framed_pentamer_count = $framed_pentamers{"${pentamer}-${frame}"};

		my $markov_prob_framed = $framed_hexamer_count / $framed_pentamer_count;
		
		## get the background freqs and prob
		my $background_hexamer_count = $background_hexamers{$hexamer};
		my $background_pentamer_count = $background_pentamers{$pentamer};

		my $background_prob = $background_hexamer_count / $background_pentamer_count;

		my $logliklihood = log($markov_prob_framed / $background_prob);

		print "$framed_hexamer\t$logliklihood\n";
	}





	## The Initialization Matrix based on framed pentamer frequencies.

	## get total pentamers
	my $total_pentamer_counts = 0;
	foreach my $background_pentamer (keys %background_pentamers) {
		
		my $background_count = $background_pentamers{$background_pentamer};
		$total_pentamer_counts += $background_count;
	}
	

	foreach my $framed_pentamer (sort keys %framed_pentamers) {
		
		my ($pentamer, $frame) = split (/-/, $framed_pentamer);

		my $frame_counts = $framed_all_pentamer_counts{$frame};
		my $framed_pentamer_counts = $framed_pentamers{$framed_pentamer};

		my $prob_framed_pentamer = $framed_pentamer_counts / $frame_counts;
		
		## now background
		my $background_pentamer_counts = $background_pentamers{$pentamer};
		my $prob_background_pentamer = $background_pentamer_counts / $total_pentamer_counts;

		my $logliklihood = log($prob_framed_pentamer / $prob_background_pentamer);

		print "$framed_pentamer\t$logliklihood\n";

	}

	return;
}

####
sub add_pseudocounts {
	
	foreach my $framed_hexamer (keys %framed_hexamers) {
		my ($hexamer, $frame) = split (/-/, $framed_hexamer);
		
		my $pentamer = substr($hexamer, 0, 5);
		
		$framed_hexamers{$framed_hexamer}++;
		$framed_pentamers{"${pentamer}-${frame}"}++;
		$framed_all_pentamer_counts{$frame}++;

		$background_hexamers{$hexamer}++;
		$background_pentamers{$pentamer}++;

	}


	return;
}

	
####
sub parse_targetCDSs {
	my ($seqFile) = @_;

	my $fasta_reader = new Fasta_reader($seqFile);
	
	while (my $seq_obj = $fasta_reader->next()) {
		
		my $accession = $seq_obj->get_accession();
		print STDERR "\r     Target: processing $accession           ";
		
		my $sequence = uc $seq_obj->get_sequence();

		my $seq_len = length($sequence);

		for (my $i = 0; $i <= $seq_len - 5; $i++) {
			my $frame = $i % 3;
			my $pentamer = substr($sequence, $i, 5);
			$framed_pentamers{"${pentamer}-${frame}"}++;
			$framed_all_pentamer_counts{$frame}++;
			
			if ($i <= $seq_len - 6) { 
				# got a hexamer
				my $hexamer = substr($sequence, $i, 6);
				$framed_hexamers{"${hexamer}-${frame}"}++;
			}
		}
	}
	
	return;
}

#### 
sub parse_background {
	my ($seqFile) = @_;

	my $fasta_reader = new Fasta_reader($seqFile);                                                                                                    
	
    while (my $seq_obj = $fasta_reader->next()) {
		
		my $accession = $seq_obj->get_accession();

		print STDERR "\r   background: processing $accession           ";
		
		my $sequence = uc $seq_obj->get_sequence();		
        my $seq_len = length($sequence);    

		## do the forward and reverse complement

		foreach my $seq ($sequence, &reverse_complement($sequence)) {
			
			for (my $i = 0; $i <= $seq_len - 5; $i++) {                                                                                                          
				
				my $pentamer = substr($seq, $i, 5);                                                                                                         
				$background_pentamers{"$pentamer"}++;                                                                                                       
				
				if ($i <= $seq_len - 6) {                                                                                                                        
					# got a hexamer                                                                                                                              
					my $hexamer = substr($sequence, $i, 6);                                                                                                      
					$background_hexamers{"$hexamer"}++;                                                                                                     
				}                                                                                                                                                
			}                                                                                                                                                    
		}                                                    
	}
	
	return;
}

