#!/usr/bin/env perl

use strict;
use warnings;

use Carp;

## minion lastal settings: lastal -q 1 -a 1 -b 1 

my $usage = "usage: $0 last.output [incl_pretty_aln]\n\n";

my $last_file = $ARGV[0] or die $usage;
my $include_pretty_align_flag = $ARGV[1] || 0;


main: {

    ## print column header
    print join("\t", "#query", "query_coords", "db_hit", "db_hit_coords", 
               "num_aligned_chars", "per_ID", "per_gap_query", "per_gap_db_hit", 
               "matches", "mismatches", "gap_query", "gap_db_hit") . "\n";
    

    open (my $fh, $last_file) or die $!;
    while (<$fh>) {
        chomp;
        if (/^\#/) { next; }
        unless (/\w/) { next; }
        
        if (/^a score=(\d+)/) {
            my $score = $1;
            
            my $db_entry = <$fh>;
            chomp $db_entry;
            
            my $query_entry = <$fh>;
            chomp $query_entry;

            &report_hit($score, $db_entry, $query_entry);
        }
        
    }

    exit(0);
}

####
sub report_hit {
    my ($score, $db_entry, $query_entry) = @_;

    my $db_struct = &parse_align($db_entry);
    my $query_struct = &parse_align($query_entry);

    my $align_stats = &examine_alignment($db_struct->{seq}, $query_struct->{seq});
    
    print join("\t", $query_struct->{acc}, "$query_struct->{lend}-$query_struct->{rend}\[$query_struct->{orient}]", 
               $db_struct->{acc}, "$db_struct->{lend}-$db_struct->{rend}\[$db_struct->{orient}]",
               
               $align_stats->{aligned}, $align_stats->{per_ID}, $align_stats->{per_gap_A}, $align_stats->{per_gap_B},
               $align_stats->{matches}, $align_stats->{mismatches}, $align_stats->{gap_A}, $align_stats->{gap_B}) . "\n";
    

    if ($include_pretty_align_flag) {
        &pretty_alignment($db_struct, $query_struct);
    }

    return;
   
}

####
sub parse_align {
    my ($line) = @_;

    
    # MAF format described here:
    # http://genome.ucsc.edu/FAQ/FAQformat.html#format5
    
    my @x = split(/\s+/, $line);
    my ($s, $acc, $lend, $len, $orient, $size, $aligned_sequence) = @x;

    my $gapless_seq = $aligned_sequence;
    $gapless_seq =~ s/\-//g;
    
    my $align_seq_len = length($gapless_seq);
    
    $lend++; # 1-based coordinate now.
    if ($orient eq '-') {
        # revcomp the coordinate
        $lend = $size - $lend + 1;
    }
    my $rend = $lend + $align_seq_len - 1;

    my $struct = { acc => $acc,
                   
                   lend => $lend,
                   rend => $rend,
                   len => $len,
                   orient => $orient,

                   seq => $aligned_sequence,
               };

    return($struct);
}


####
sub examine_alignment {
    my ($A_seq, $B_seq) = @_;

    if (length($A_seq) != length($B_seq)) {
        confess "Error, align seq lengths aren't equal";
    }
    
    my $matches = 0;
    my $mismatches = 0;
    my $gap_A = 0;
    my $gap_B = 0;

    my @A_chars = split(//, uc $A_seq);
    my @B_chars = split(//, uc $B_seq);

    for (my $i = 0; $i < $#A_chars; $i++) {
        
        my $A_char = $A_chars[$i];
        my $B_char = $B_chars[$i];
        
        if ($A_char eq '-') {
            $gap_A++;
        }
        elsif ($B_char eq '-') {
            $gap_B++;
        }
        elsif ($A_char eq $B_char) {
            $matches++;
        }
        elsif ($A_char ne $B_char) {
            $mismatches++;
        }
        else {
            die "Error, shouldn't get here...";
        }

    }
    
    my $aligned_chars = $matches + $mismatches + $gap_A + $gap_B;
    
    my $struct = {  matches => $matches,
                    mismatches => $mismatches,
                    
                    aligned => $aligned_chars,
                    
                    gap_A => $gap_A,
                    gap_B => $gap_B,

                    per_ID => sprintf("%.2f", $matches/$aligned_chars * 100),
                    
                    per_gap_A => sprintf("%.2f", $gap_A / $aligned_chars * 100),
                    per_gap_B => sprintf("%.2f", $gap_B / $aligned_chars * 100),
                };

    return($struct);
}
                    
        
        
####
sub pretty_alignment {
    my ($db_struct, $query_struct) = @_;

    my $db_acc = $db_struct->{acc};
    my $query_acc = $query_struct->{acc};

    my $max_acc_len = length($db_acc);
    if (length($query_acc) > $max_acc_len) {
        $max_acc_len = length($query_acc);
    }
    
    my $align_length = 60;

    my $db_aln_seq = $db_struct->{seq};
    my $query_aln_seq = $query_struct->{seq};



    my $idx = 0;
    while ($idx < length($db_aln_seq)) {
        
        my $db_seq_region = uc substr($db_aln_seq, $idx, $align_length);
        my $query_seq_region = uc substr($query_aln_seq, $idx, $align_length);

        my @db_chars = split(//, $db_seq_region);
        my @query_chars = split(//, $query_seq_region);
        my @match_chars = ();

        for (my $i = 0; $i <= $#db_chars; $i++) {
            if ($db_chars[$i] eq $query_chars[$i]) {
                push (@match_chars, '|');
            }
            else {
                push (@match_chars, ' ');
            }
        }

        print sprintf("%-${max_acc_len}s  ", $db_acc) . join("", @db_chars) . "\n";
        print sprintf("%-${max_acc_len}s  ", "") . join("", @match_chars) . "\n";
        print sprintf("%-${max_acc_len}s  ", $query_acc) . join("", @query_chars) . "\n";
        print "\n";
        
        $idx += $align_length;
    
    }

    return;

}
