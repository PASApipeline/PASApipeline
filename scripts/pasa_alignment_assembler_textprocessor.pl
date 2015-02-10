#!/usr/bin/env perl

use FindBin;
use lib ($FindBin::Bin);
use Pasa_init;
use lib ("../PerlLib");
use strict;
use CDNA::PASA_alignment_assembler;
use CDNA::CDNA_alignment;
use Getopt::Std;

use vars qw ($opt_h $opt_D $opt_p $opt_d $DEBUG $opt_S $opt_M $opt_s);

&getopts ('hD:dp:S:M:s');


$|=1;
our $SEE = $opt_s;

open (STDERR, "&>STDOUT");

my $usage =  <<_EOH_;

usage: $0 < alignment.textfile

The alignment.textfile should have the following format:

cdna_acc,orient,full-length_flag,segment_coords,...

If the cDNA is full-length, full-length_flag = 1, else 0

The strand should be '+|-|?'  
Use '?' in cases where you have single exon alignments of ambiguous transcribed orientation.

ie. 

// cluster: 130
gi|20259098|gb|AY091326.1|,-,1,30806-31070,31156-31467,31670-31809,31970-32033,32122-32553,32699-32836,33131-34534
gi|18377677|gb|AY074292.1|,-,1,30611-31070,31156-31467,31670-31809,31970-32033,32122-32553,32699-32836,33131-34545,34644-35373,35682-36092
AI996609,-,0,30641-31070,31156-31247
AV439477,+,0,30742-31027
AV441677,-,0,31422-31467,31670-31809,31970-32033,32122-32416
AV522465,-,0,30613-31070,31156-31223
AV522962,-,0,30684-31070,31156-31381


// cluster ... etc, etc...

Namely,
accession,orientation,coordinates

_EOH_

    ;

if ($opt_h) { die $usage;}


$/ = "\n//";

while (my $input = <STDIN>) {
    print "######################################################\n";
    print $input . "\n";
    my $assembler = new CDNA::PASA_alignment_assembler();
    
    my @datalines = split (/\n/, $input);
    my $header = "";
    if ($datalines[0] =~ /\/\//) {
        $header = shift @datalines; #lose the // line
    }
    chomp $header;
    if ($datalines[$#datalines] =~ /\/\//) { pop @datalines;} #rid the last //-containing line.
    my @alignments;
    
    my %seen;
   
    foreach my $dataline (@datalines) {
        unless ($dataline =~ /,/) { next;}
        my ($acc, $strand, @coordsets) = split (/,/, $dataline);
        
        if ($seen{$acc}) {
            die "Error, cannot report alignments for acc($acc) multiple times in a single input to pasa.\n";
        }
        $seen{$acc} = 1;
        
        my @segments;
        my $cdna_length = 0;
        foreach my $coordset (@coordsets) {
            my ($lend, $rend) = split (/-/, $coordset);
            my $segment = new CDNA::Alignment_segment($lend, $rend);
            push (@segments, $segment);
            $cdna_length += abs ($rend - $lend) + 1;
        }
        my $alignment = new CDNA::CDNA_alignment($cdna_length, \@segments);
        if ($strand =~ /^[\+\-]$/) {
            $alignment->force_spliced_validation($strand);
        } 
        elsif ($strand eq "?") {
            $alignment->set_spliced_orientation($strand);
        }
        
        
        $alignment->set_acc($acc);
        push (@alignments, $alignment);
    }
    if (@alignments) {
        print "HEADER: $header\n" if $header;
        $assembler->assemble_alignments(@alignments);
        # set orientation for display purposes.
        my @assemblies = $assembler->get_assemblies();
        foreach my $assembly (@assemblies) {
            if ((my $orient = $assembly->get_spliced_orientation()) ne '?') {
                $assembly->force_spliced_validation($orient);
            }
            $assembly->remap_cdna_segment_coords();
        }
        
        print $assembler->toAlignIllustration(60);
        print "\n\n\n";
    }
    my $x=0;
    foreach my $assembly ($assembler->get_assemblies()) {
        $x++;
        print "Assembly($x): " .  $assembly->toToken() . "\n";
    }
}

exit(0);


