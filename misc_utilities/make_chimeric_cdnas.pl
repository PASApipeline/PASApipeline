#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;

my $usage = "usage: $0 cdna.fasta num_chimeras\n\n";

my $cdna_fasta = $ARGV[0] or die $usage;
my $num_chimeras = $ARGV[1] or die $usage;

main: {
    
    my $fasta_reader = new Fasta_reader($cdna_fasta);
    
    my %trans_to_info;

    while (my $seq_obj = $fasta_reader->next()) {


        my $header = $seq_obj->get_header();

        my $sequence = $seq_obj->get_sequence();

        #print "header: $header\nseq: $sequence\n\n";
        
        my @exon_coords = &get_exon_coords($sequence);
        
        if (0) {
            foreach my $coordset (@exon_coords) {
                print join("-", @$coordset) . "\t";
                my $subseq = substr($sequence, $coordset->[0], $coordset->[1] - $coordset->[0] + 1);
                print "$subseq\n";
            }
            print "\n\n";
        }
    

        unless (scalar @exon_coords > 1) { next; } ## only multi-exonic ones for now.
    
        my ($trans, $gene) = split(/\s+/, $header);
        
        my $struct = { gene => $gene,
                       trans => $trans,
                       
                       sequence => $sequence,
                       
                       exon_coords => [@exon_coords],
        };

        $trans_to_info{$trans} = $struct;


    }


    ###
    # Randomly make chimeras
    ###

    my @trans = keys %trans_to_info;
    my $num_trans = scalar(@trans);

    my $chimera_counter = 0;
    while ($chimera_counter < $num_chimeras) {

        my $trans_A_struct = $trans_to_info{ $trans[ int(rand($num_trans)) ] };
        my $trans_B_struct = $trans_to_info{ $trans[ int(rand($num_trans)) ] };

        if ($trans_A_struct->{gene} eq $trans_B_struct->{gene}) {
            print STDERR "... skipping $trans_A_struct->{gene}, same gene in both selections.\n";
            next;
            
        }

        my $trans_A = $trans_A_struct->{trans};
        my $sequence_A = $trans_A_struct->{sequence};        
        my @trans_A_coords = @{$trans_A_struct->{exon_coords}};
        my $num_exons_A = scalar(@trans_A_coords);

        my $trans_B = $trans_B_struct->{trans};
        my $sequence_B = $trans_B_struct->{sequence};
        my @trans_B_coords = @{$trans_B_struct->{exon_coords}};
        my $num_exons_B = scalar(@trans_B_coords);

        my $clip_pt_A = $trans_A_coords[ int(rand($num_exons_A-1)) ]->[1]; # all but the last exon are fair game
        my $clip_pt_B = $trans_B_coords[ int(rand($num_exons_B-1)) + 1]->[0]; # all but the first exon are fair game
        
        
        my $chimera_A_seq = substr($sequence_A, 0, $clip_pt_A+1);
        #print "SeqA:\n$sequence_A\nptA:\n$chimera_A_seq\n\n";
        
        my $chimera_B_seq = substr($sequence_B, $clip_pt_B);
        #print "SeqB:\n$sequence_B\nptB:\n$chimera_B_seq\n\n";
        
        print ">${trans_A}--${trans_B}\n" . uc($chimera_A_seq) . lc($chimera_B_seq) . "\n";
        
        
        $chimera_counter++;
    }


    exit(0);

}



####
sub get_exon_coords {
    my ($sequence) = @_;

    my $seq_length = length($sequence);

    my @coords;
    while ($sequence =~ /([A-Z]+)/g) {
        my $beg = $-[0];
        my $end = $+[0];
        $end--;
    
        if (! @coords) {
            # first one
            if ($beg != 0) {
                push (@coords, [0, $beg-1]);
            }
        }
        else {
            # not first one
            my $last_pos = $coords[$#coords]->[1];
            my $internal_beg = $last_pos + 1;
            my $internal_end = $beg - 1;
            push (@coords, [$internal_beg, $internal_end]);
        }
        push (@coords, [$beg, $end]);
        
    }
    
    my $last_pos = $coords[$#coords]->[1];
    if ($last_pos != $seq_length - 1) {
        push (@coords, [$last_pos + 1, $seq_length - 1]);
    }

    return(@coords);
    
}


