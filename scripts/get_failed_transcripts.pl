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
use Getopt::Std;
use Fasta_reader;

use vars qw ($opt_c $opt_t);

&getopts ('c:t:');

$|=1;
our $SEE = 0;

open (STDERR, "&>STDOUT");

my $usage =  <<_EOH_;

############################# Options ###############################
#
# -c * <filename>               configuration file for align-assembly
# -t transcripts <filename>     
#
###################### Process Args and Options #####################

_EOH_

    ;


my $configfile = $opt_c or die $usage;
my $transcript_db = $opt_t or die $usage;

## Read configuration file.
my %config = &readConfig($configfile);


my $mysql_db = $config{MYSQLDB} or die "Error, couldn't extract mysql_db name from config file " . cwd() . "/$configfile\n";
my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $user = &Pasa_conf::getParam("MYSQL_RW_USER");
my $password = &Pasa_conf::getParam("MYSQL_RW_PASSWORD");

## Create the database
my $admin_db = ""; # just forcing a connection.
my ($dbproc) = &connect_to_db($mysql_server,$admin_db,$user,$password);

my $query = "use $mysql_db";
&RunMod($dbproc, $query);


$query = "select ci.cdna_acc from cdna_info ci where ci.is_assembly = 0 and not exists (select 1 from align_link al where al.cdna_info_id = ci.id and al.validate = 1)";
my @results = &do_sql_2D($dbproc, $query);

my %failed_accs;
my $num_failed = 0;
foreach my $result (@results) {
    #print $result->[0] . "\n";
    
    my $acc = $result->[0];
    $failed_accs{$acc} = 1;
    $num_failed++;
}

print STDERR "-identified $num_failed failed alignments.\n";

my $fasta_reader = new Fasta_reader($transcript_db);
while (my $seq_obj = $fasta_reader->next()) {
    
    my $acc = $seq_obj->get_accession();
    if ($failed_accs{$acc}) {
        
        my $entry = $seq_obj->get_FASTA_format();
        print $entry;
    }
}


$dbproc->disconnect;

exit(0);
