#!/usr/bin/env perl

use strict;
use FindBin;
use lib ($FindBin::Bin);
use Pasa_init;
use CDNA::Overlap_assembler;

our $SEE = 0;

my %assemblers;

if (-d "clusters") {
    die "Please remove the existing \'clusters\' directory before continuing.";
} else {
    mkdir "clusters";
    chmod (0775, "clusters");
}

while (<STDIN>) {
    my @x = split (/\t/);
    unless ($x[0]=~ /^\d+$/) { print STDERR "skipping $_" if $SEE; next;}
    my ($cdna_acc, $chromo_acc, $cdna_start, $cdna_stop) = ($x[9], $x[13], $x[15], $x[16]);
    my $assembler;
    unless ($assembler = $assemblers{$chromo_acc}) {
	print "Creating assembler for $chromo_acc.\n" if $SEE;
	$assembler = new CDNA::Overlap_assembler();
	$assemblers{$chromo_acc} = $assembler;
    }
    $assembler->add_cDNA($cdna_acc, $cdna_start, $cdna_stop);
}


print STDERR "Done building data structure.\n" if $SEE;

## build clusters for each chromosome:
foreach my $chromo (sort keys %assemblers) {
    ## order and assign index to nodes.
    my $assembler = $assemblers{$chromo};
    
    ## peform single-linkage-clustering to group together all overlapping cDNAs:
    print "\n\n// Processing $chromo clusters.\n" if $SEE;
    my @clusters = $assembler->build_clusters();
    open (FILE, ">clusters/$chromo.clusters") or die "Can't open file here.\n";
    foreach my $set_ref (@clusters) {
	my $line = join ("\t", @$set_ref);
	print FILE "$line\n";
    }
    close FILE;
}

