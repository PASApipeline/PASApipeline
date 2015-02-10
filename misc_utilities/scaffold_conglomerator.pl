#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");         
use Fasta_reader;

my $usage = "usage: $0 genomeFastaFile numConglomerates prefix\n\n";

my $genome_file = $ARGV[0] or die $usage;
my $numConglomerates = $ARGV[1] or die $usage;
my $prefix = $ARGV[2] or die $usage;

main: {
    print STDERR "-computing total genome size.\n";
    my $totalSeqLength = &sum_seq_lengths($genome_file);
    print STDERR "\n\nTotal sequence length of $genome_file = $totalSeqLength\n";
    
    ## open output files:
    open (my $seqFH, ">$prefix.conglomerates.fasta") or die $!;
    open (my $scaffFH, ">$prefix.conglomerates.scaffolding") or die $!;
    
    my $max_mol_length = int($totalSeqLength / $numConglomerates + 0.5);
    
    my $conglomerate_seq = "";
    my $start_pt = 0;
    my $mol_count = 1;

    my $fasta_reader = new Fasta_reader($genome_file);
    while (my $seqObj = $fasta_reader->next()) {
        my $genome_seq = $seqObj->get_sequence();
        my $acc = $seqObj->get_accession();

        my $len = length($genome_seq);
        
        print STDERR "-processing $acc, length: $len, currently: $mol_count / $numConglomerates\n";
    
        ## append to conglomerate
        $conglomerate_seq .= $genome_seq;
        my $begin = $start_pt + 1;
        my $end = $start_pt + $len;
        print $scaffFH "${prefix}_p$mol_count\t$acc\t$begin\t$end\tlen: $len\n";
        $start_pt = $end;
    
        if ($end > $max_mol_length) { # not really a max then, is it.  ;)
            ## dump it, start next mol:
            $conglomerate_seq =~ s/(\S{60})/$1\n/g;
            chomp $conglomerate_seq;
            print $seqFH ">${prefix}_p${mol_count}\n$conglomerate_seq\n";
            
            # next mol
            $mol_count++;
            
            # reinit
            $start_pt = 0;
            $conglomerate_seq = "";
        }
    }
    
    if ($conglomerate_seq) {
        # write last conglomerate:
        $conglomerate_seq =~ s/(\S{60})/$1\n/g;
        chomp $conglomerate_seq;
        print $seqFH ">${prefix}_p${mol_count}\n";
    }

    close $seqFH;
    close $scaffFH;
    
    exit(0);
}


####
sub sum_seq_lengths {
    my ($genome_file) = @_;
    
    my $sum_length = 0;

    my $fasta_reader = new Fasta_reader($genome_file);
    while (my $seqObj = $fasta_reader->next()) {
        my $genome_sequence = $seqObj->get_sequence();
        my $acc = $seqObj->get_accession();
        my $len = length($genome_sequence);

        $sum_length += $len;

        print STDERR "\rread $acc, len=$len, sum=$sum_length                     ";
        
    }

    return ($sum_length);
}
