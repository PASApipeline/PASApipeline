#!/usr/bin/env perl

use strict;
use warnings;

my $seg_per_id = $ARGV[0] || 0;


my $transcript_acc = "";
my $chain_counter = 0;
my $segment_counter = 0;

while (<STDIN>) {
    if (/^>(\S+)/) {
        $transcript_acc = $1;
        
    }
	
	elsif (/Alignment for path \d+:/) {
		$chain_counter++;
        $segment_counter = 0; #reinit
	}
	
    elsif (/\s*([\+\-])([^\:]+):(\d+)-(\d+)\s+\((\d+)-(\d+)\)\s+(\d[^\%]+\%)/) {
        #print " $transcript_acc\t$_";
        my $orient = $1;
        my $genome_acc = $2;
        my $genome_end5 = $3;
        my $genome_end3 = $4;
        my $transcript_end5 = $5;
        my $transcript_end3 = $6;
        my $percent_identity = $7;
        
        $percent_identity =~ s/\%//;
        
        $segment_counter++;
        my @x;
        $#x = 14; #prealloc
        # init list
        foreach my $ele (@x) {
            $ele = "";
        }
        
        $x[0] = $genome_acc;
        $x[3] = "gmap";
        $x[5] = $transcript_acc;
        $x[6] = $genome_end5;
        $x[7] = $genome_end3;
        $x[8] = $transcript_end5;
        $x[9] = $transcript_end3;
        $x[10] = $percent_identity;
        $x[13] = $chain_counter;
        $x[14] = $segment_counter;
        
        my $btab_line = join ("\t", @x);
        print "$btab_line\n" if ($percent_identity >= $seg_per_id);
        
    }
}

exit(0);

