#!/usr/bin/env perl

use FindBin;
use lib ($FindBin::Bin);
use Pasa_init;
use Mysql_connect;
use strict;
use DBI;
use Getopt::Std;
use Ath1_cdnas;
use CDNA::CDNA_alignment;
use CDNA::Alternative_splice_comparer;
use CDNA::PASA_alignment_assembler;
use Carp;
use Data::Dumper;
use CdbTools;
use Nuc_translator;

use vars qw ($opt_M $opt_v $opt_p $opt_G $opt_d $opt_h $opt_g);
open (STDERR, "&>STDOUT");
&getopts ('M:p:dhvg:');
my $usage =  <<_EOH_;

Script loads the alignment textual representation for the pasa assemblies.

############################# Options ###############################
# -M Mysql database/server ie. ("ath1_cdnas:haasbox")
# -p passwordinfo  (contains "username:password")
# -g gff3 file
# -d Debug
# 
# -h print this option menu and quit
# -v verbose
###################### Process Args and Options #####################

_EOH_

    ;

my $SEE = $opt_v;
our $DB_SEE = $opt_v;

if ($opt_h) {die $usage;}
my $MYSQLstring = $opt_M or die $usage;
my ($MYSQLdb, $MYSQLserver) = split (/:/, $MYSQLstring); 
my $passwordinfo = $opt_p or die $usage;
our $DEBUG = $opt_d;
my $gff_file = $opt_g or die $usage;

my ($user, $password) = split (/:/, $passwordinfo);

my %cdna_info;

my $dbproc;

main: {
     
    $dbproc = &connect_to_db($MYSQLserver, $MYSQLdb, $user, $password);
    
    
    ## get cdna_length info
    &get_cdna_info();
    
    print STDERR "-parsing and loading alignments.\n";
    ## parse data from the gff file:
    &parse_GFF3_file_load_alignments();
    
    print STDERR "\nFinished.\n\n";
    
    $dbproc->disconnect;
    exit(0);
    
}




####
sub get_cdna_info {
    my $query = qq { select cdna_acc, id, length from cdna_info };
    my @results = &do_sql_2D($dbproc, $query);
    foreach my $result (@results) {
        my ($cdna_acc, $id, $length) = @$result;
        $cdna_info{$cdna_acc} = { id => $id,
                                  length => $length};
    }
}



####
sub parse_GFF3_file_load_alignments {

    my @match_list;
    my $current_match = "";
        
    my %asmbl_to_acc;

    my $counter = 0;
    open (my $fh, $gff_file) or die "Error, cannot open $gff_file ";
    while (<$fh>) {
        unless (/\w/) { next; }
		#print;
        chomp;
        my @x = split (/\t/);
        my ($asmbl_id, $end5, $end3, $per_id, $orient, $match_info) = ($x[0], $x[3], $x[4], $x[5], $x[6], $x[8]);
        
        $match_info =~ /ID=([^\s;]+);?/;
        my $match_id = $1 or die "Error, couldn't parse match from $match_info, line: $_\n";
        
        $match_info =~ /Target=([^\s;]+) (\d+) (\d+);?/;
        my $acc = $1 or die "Error, couldn't parse target from $match_info, line: $_\n";
        my $match_lend = $2;
        my $match_rend = $3;

        $acc = $match_id; ## NOTE, resetting this to include the path number (align_acc instead of cdna_acc)


        ## track asmbl,acc info so we can track this in the database.
        $asmbl_to_acc{$asmbl_id}->{$acc} = 1;
        
        unless ($match_rend =~ /^\d+$/ && $match_lend =~ /^\d+$/) {
            die "error parsing match coordinates lend[$match_lend] rend[$match_rend]from last field: $match_info, line: $_\n";
        }
        
        if ($match_id ne $current_match) {
            if (@match_list) {
                $counter++;
                # if ($counter > 300) { last;}
                print STDERR "\rloading $counter   ";
                &process_match_list (@match_list);
            }
            $current_match = $match_id;
            @match_list = ();
        }
        
        
        push (@match_list, { end5 => $end5,
                             end3 => $end3,
                             orient => $orient,
                             match_lend => $match_lend,
                             match_rend => $match_rend,
                             per_id => $per_id,
                             match_id => $match_id,
                             acc => $acc,
                         }
              );
        
        
        
    }
    close $fh;
    
    ## get last one.
    if (@match_list) {
        &process_match_list (@match_list);
    }
    
    
    &store_cluster_links(\%asmbl_to_acc);
        
}


####
sub store_cluster_links {
    my $asmbl_to_acc_href = shift;

    print "\n\nStoring cluster links.\n\n";
    
    foreach my $asmbl_id (keys %$asmbl_to_acc_href) {
        my @accs = keys %{$asmbl_to_acc_href->{$asmbl_id}};

        ## install row in clusters table

        my $query = qq { insert clusters (annotdb_asmbl_id) values (?) };
        &RunMod($dbproc, $query, $asmbl_id);
        
        my $cluster_id = &Mysql_connect::get_last_insert_id($dbproc);

        ## set the cluster link values
        $query = qq { update align_link set cluster_id = $cluster_id where align_acc = ? };
        my $sth = $dbproc->{dbh}->prepare($query);
        foreach my $acc (@accs) {
                        
            my $info = $cdna_info{$acc};

            
            $sth->execute($acc);
            
        }
        
        $sth->finish;
    }
}





#### 
sub process_match_list {
    my @matches = @_;
    
    my $acc = $matches[0]->{acc};
    my $orient = $matches[0]->{orient};
    
    my @segments;
    
    my $alignment_length = 0;
    my $sum_per_id = 0;
    
    foreach my $match (@matches) {
        my ($end5, $end3, $orient, $match_lend, $match_rend, $per_id) = ($match->{end5},
                                                                         $match->{end3},
                                                                         $match->{orient},
                                                                         $match->{match_lend},
                                                                         $match->{match_rend},
                                                                         $match->{per_id},
                                                                         );
        
        if ($orient eq '-') {
            ($end5, $end3) = ($end3, $end5);
        }
    
        $alignment_length += $match_rend - $match_lend + 1;
        $sum_per_id += $per_id;
    
        
        my $segment = new CDNA::Alignment_segment ($end5, $end3, $match_lend, $match_rend, $per_id);
        push (@segments, $segment);
    }
    
    my $real_trans_acc = $acc;
    $real_trans_acc =~ s/\.path\d+$//;
    
    my $cdna_info = $cdna_info{$real_trans_acc} or die "Error, no known cdna_length for accession($real_trans_acc) ";
    my $cdna_id = $cdna_info->{id};
    my $cdna_length = $cdna_info->{length};
    
    my $alignment = new CDNA::CDNA_alignment ($cdna_length, \@segments);
    
    $alignment->set_acc($acc);
    $alignment->set_cdna_id($cdna_id);
    $alignment->set_spliced_orientation($orient);
    $alignment->set_orientation($orient);
    
    my $align_id = &load_alignment($alignment);
    
    my $avg_per_id = sprintf ("%.2f", $sum_per_id / $alignment_length * 100);
    
    my $query = qq { update align_link set avg_per_id = ? where align_id = $align_id };

    &RunMod($dbproc, $query, $avg_per_id);
    


}
 

####
sub load_alignment {
    my ($alignment) = @_;
    
    #print "Loading: " . $alignment->get_acc() . " " . $alignment->toToken() . "\n";
    
    
    my $align_id = &Ath1_cdnas::load_CDNA_alignment_obj($dbproc, $alignment, "gmap", undef);
    
    return ($align_id);
    
}

