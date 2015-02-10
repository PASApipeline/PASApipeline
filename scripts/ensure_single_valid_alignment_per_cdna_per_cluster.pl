#!/usr/bin/env perl

use FindBin;
use lib ($FindBin::Bin);
use Pasa_init;
use Pasa_conf;
use Mysql_connect;
use strict;
use ConfigFileReader;
use DBI;
use Cwd;
use Ath1_cdnas;
use Getopt::Long qw (:config no_ignore_case bundling pass_through);


our $SEE = 0;


my $usage = <<__EOUSAGE__;

###############################################################################
#
# -M <string>     mysql database name
#
#
###############################################################################


__EOUSAGE__

    ;

my $mysql_db;

&GetOptions("M=s" => \$mysql_db);

unless ($mysql_db) {
    die $usage;
}


my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $user = &Pasa_conf::getParam("MYSQL_RW_USER");
my $password = &Pasa_conf::getParam("MYSQL_RW_PASSWORD");

my ($dbproc) = &connect_to_db($mysql_server,$mysql_db,$user,$password);

my $query = "use $mysql_db";
&RunMod($dbproc, $query);

my $query = "select ci.cdna_acc, a.cluster_id, a.cdna_info_id, a.align_id, a.align_acc, a.prog, a.score "
    . " from cdna_info ci, align_link a "
    . " where ci.id = a.cdna_info_id "
    . " and a.validate = 1 "
    . " and ci.is_assembly = 0 "
    . " order by cluster_id, cdna_info_id";
    
my @results = &do_sql_2D($dbproc, $query);

my %cluster_n_cdna_info;

foreach my $result (@results) {
    
    my ($cdna_acc, $cluster_id, $cdna_info_id, $align_id, $align_acc, $prog, $score) = @$result;
    
    my $token = join("$;", $cluster_id, $cdna_info_id);
    push (@{$cluster_n_cdna_info{$token}}, { 
        cdna_acc => $cdna_acc,
        align_acc => $align_acc,
        align_id => $align_id,
        prog => $prog,
        score => $score,
        cluster_id => $cluster_id,
    });
    
}

foreach my $collection_aref (values %cluster_n_cdna_info) {
    
    my @alignments = @{$collection_aref};
    
    if (scalar @alignments > 1) {
        
        @alignments = reverse sort {$a->{score}<=>$b->{score}} @alignments;
        
        my $best_alignment = shift @alignments;
        
        print "KEEPING:\t" . join("\t", 
                                  $best_alignment->{cluster_id},
                                  $best_alignment->{cdna_acc},
                                  $best_alignment->{align_acc},
                                  $best_alignment->{align_id},
                                  $best_alignment->{prog},
                                  $best_alignment->{score}) . "\n";
        
        ## invalidate the others
        foreach my $remaining_alignment (@alignments) {
            
            print "INVALIDATING:\t" . join("\t", 
                                           $remaining_alignment->{cluster_id},
                                           $remaining_alignment->{cdna_acc},
                                           $remaining_alignment->{align_acc},
                                           $remaining_alignment->{align_id},
                                           $remaining_alignment->{prog},
                                           $remaining_alignment->{score}) . "\n";
            
            my $query = "update align_link set comment = 'invalidating in preference for diff alignment at location', validate = 0 where align_id = " . $remaining_alignment->{align_id};

            &RunMod($dbproc, $query);
        
        }
        print "\n"; # spacer
        
    }
    
}



$dbproc->disconnect;

exit(0);
