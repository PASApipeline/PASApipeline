#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my $usage = "usage: $0 top_hits_btab_1 top_hits_btab_2 (btab|m2fmt)\n\n";

my $top_hits_1_file = $ARGV[0] or die $usage;
my $top_hits_2_file = $ARGV[1] or die $usage;
my $format = $ARGV[2] or die $usage;

 main: {
     
     my %top_hits_1 = &get_top_hits($top_hits_1_file);
     my %top_hits_2 = &get_top_hits($top_hits_2_file);
     
     foreach my $acc_1 (keys %top_hits_1) {
         my $best_hit_acc = $top_hits_1{$acc_1}->{acc};
         my $best_hit_e_value = $top_hits_1{$acc_1}->{e_value};
		 my $best_hit_per_id = $top_hits_1{$acc_1}->{per_id};

         my $best_hit_2_ref = $top_hits_2{$best_hit_acc};
         
         if (ref $best_hit_2_ref) {
             my $best_hit_acc2 = $best_hit_2_ref->{acc};
             my $best_hit_evalue2 = $best_hit_2_ref->{e_value};
             my $best_hit_per_id2 = $best_hit_2_ref->{per_id};

             #print "$acc_1\t$best_hit_acc2\n";
             
             if ($acc_1 eq $best_hit_acc2) {
                 print "ortholog?\t$acc_1\t$best_hit_acc\tAvsB: $best_hit_e_value\t$best_hit_per_id\tBvsA: $best_hit_evalue2\t$best_hit_per_id2\n";
             }
         }
     }
          
 }


exit(0);


####
sub get_top_hits {
    my ($file) = @_;
    
    

    my %top_hits;
    open (my $fh, $file) or die $!;
    while (<$fh>) {
        my @x = split (/\t/);
		
		my ($query_acc, $db_acc, $per_id, $e_value);

		if ($format eq 'btab') {
			$per_id = $x[10];
			$query_acc = $x[0];
			$db_acc = $x[5];
			$e_value = $x[19];
		}
		elsif ($format eq 'm2fmt') {
			$per_id = $x[10];
			$query_acc = $x[0];
			$db_acc = $x[1];
			$e_value = $x[2];
		}
		else {
			die "Error, don't understand format $format";
		}

		if (defined($per_id) && $per_id > 1) {
						
            if ( (exists $top_hits{$query_acc}) && $top_hits{$query_acc}->{e_value} < $e_value) { 
                # already have a better match
                next; 
            }
            
            $top_hits{$query_acc} = { acc => $db_acc,
                                      e_value => $e_value,
									  per_id => $per_id,
								  };
        }
    }
    close $fh;
    
    return (%top_hits);
}

