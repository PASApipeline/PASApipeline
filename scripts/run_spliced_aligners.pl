#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use FindBin;
use threads;

use lib ("$FindBin::Bin/../PerlLib");
use Thread_helper;
use Process_cmd;

my $usage =  <<__EOUSAGE__;

########################################################################
#
# --aligners <string>        comma-delimited list of aligners to execute
#                            Currently supported: gmap, blat
#
# --transcripts <string>     transcript database in fasta format
#
# --genome <string>          genome database in fasta format
#
# -I <int>                   maximum intron length
#
# --CPU <int>                max threads for each tool (default: 1)
#
# -N <int>                   number of top hits (default: 1)
#
#########################################################################


__EOUSAGE__

    ;

my $help_flag;
my $aligners;
my $transcripts_db;
my $genome_db;
my $max_intron_length;
my $CPU = 1;
my $num_top_hits = 1;

my %SUPPORTED_ALIGNERS = map { + $_ => 1 } qw(blat gmap);

&GetOptions ( 'h' => \$help_flag,
              'aligners=s' => \$aligners,
              'transcripts=s' => \$transcripts_db,
              'genome=s' => \$genome_db,
              'I=i' => \$max_intron_length,
              'CPU=i' => \$CPU,
              'N=i' => \$num_top_hits,
              
    );


unless ($aligners && $transcripts_db && $genome_db) {
    die $usage;
}


main: {

    my $thread_helper = new Thread_helper($CPU);
    
    if ($aligners =~ /gmap/i) {

        my $thread = threads->create('run_gmap');
        $thread_helper->add_thread($thread);
    
    }
    if ($aligners =~ /blat/i) {
        $thread_helper->wait_for_open_thread();
        my $thread = threads->create('run_blat');
        $thread_helper->add_thread($thread);
    }
    

    $thread_helper->wait_for_all_threads_to_complete();
    
    my @failed_threads = $thread_helper->get_failed_threads();
    if (@failed_threads) {
        die "Error, " . scalar(@failed_threads) . " threads failed.\n";
        exit(1);
    }
    else {
        print STDERR "processes completed successfully.\n";
        exit(0);
    }
       
}

####
sub run_gmap {
    
    my $cmd = "$FindBin::Bin/process_GMAP_alignments_gff3_chimeras_ok.pl --genome $genome_db "
        . " --transcripts $transcripts_db --CPU $CPU -N $num_top_hits ";
    if ($max_intron_length) {
        $cmd .= " -I $max_intron_length ";
    }

    $cmd .= " > gmap.spliced_alignments.gff3";

    my $checkpoint = "gmap.spliced_alignments.gff3.completed";
    
    unless (-e $checkpoint) {
        &process_cmd($cmd);
        &process_cmd("touch $checkpoint");
    }
    

    return;
}

####
sub run_blat {

    my $cmd = "$FindBin::Bin/process_BLAT_alignments.pl -g $genome_db "
        . " -t $transcripts_db -I $max_intron_length -o blat.spliced_alignments -N $num_top_hits --CPU $CPU ";

    ## checkpoints built-in to the blat runner already
    
    &process_cmd($cmd);
    
    return;
    
}
    
