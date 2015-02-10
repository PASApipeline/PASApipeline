#!/usr/bin/env perl

use FindBin;
use lib ($FindBin::Bin);
use Pasa_init;
use Pasa_conf;
use strict;
use Carp;
use File::Basename;
use Fasta_reader;
umask(0000);

my $usage = "usage: $0 genome_db transcript_db\n\nThe genome and transcript database must exist in your current directory.";

my $genome_db = $ARGV[0] or die $usage;
my $transcript_db = $ARGV[1] or die $usage;


## Set up gmap searchable database
unless (-s $genome_db && -s $transcript_db) {
    die "Error, $genome_db genome and $transcript_db transcript database must reside in the current directory.\n";
}

unless (-d "gmap_db_dir") {
    mkdir ("gmap_db_dir") or die "Error, cannot mkdir gmap_db_dir";
}

my $cmd = "gmap_setup -S -D gmap_db_dir -d $genome_db $genome_db";
&process_cmd($cmd);

$cmd = "make -f Makefile.$genome_db coords";
&process_cmd($cmd);

$cmd = "make -f Makefile.$genome_db gmapdb";
&process_cmd($cmd);

$cmd = "make -f Makefile.$genome_db install";
&process_cmd($cmd);


## Search the database one query (transcript) at a time:

my $fasta_reader = new Fasta_reader($transcript_db);

my $failures_log = "gmap.$$.failures.fasta";
open (my $failures_fh, ">$failures_log") or die $!;

my $count=0;

my $got_failure = 0;

while (my $seqobj = $fasta_reader->next()) {

	$count++;
    my $acc = $seqobj->get_accession();
    print STDERR "processing ($count) $acc\n";

    my $fasta_entry = $seqobj->get_FASTA_format();

    my $tmpfile = ".tmpfile.$$";

    open (my $fh, ">$tmpfile") or die $!;
    print $fh $fasta_entry;
    close $fh;

    my $cmd = "gmap  -D gmap_db_dir -d $genome_db -n 1 -S $tmpfile ";
    
    system $cmd;

    if ($?) {
        print STDERR "CMD $cmd died with ret $?";
        print $failures_fh $fasta_entry . "\n";
        $got_failure = 1;
    }

}

unless ($got_failure) {
    unlink $failures_log;
}


exit($got_failure);


####
sub process_cmd {
    my $cmd = shift;
    print "CMD: $cmd\n";
    my $ret = system $cmd;
    if ($ret) {
        confess "Error, cmd ($cmd) died with ret ($ret)";
    }
    else {
        return ($ret);
    }
}


