#!/usr/bin/env perl

use FindBin;
use lib ($FindBin::Bin);
use Pasa_init;
use strict;
use warnings;
use threads;
use POSIX qw(ceil);
use File::Basename;

use Fasta_reader;
use Pipeliner;

use Getopt::Long qw (:config no_ignore_case bundling);
use vars qw ($DEBUG $opt_h $opt_g $opt_t $opt_c $opt_o $opt_B $opt_I);

my $CPU = 1;
my $num_top_hits = 1;
my $KEEP_PSLX = 0;
my $output_filename = "";

&GetOptions( 'g=s' => \$opt_g,
             'd' => \$DEBUG,
             'h' => \$opt_h,
             'c=s' => \$opt_c,
             't=s' => \$opt_t,
			 'I=i' => \$opt_I,
             'N=i' => \$num_top_hits,
             'CPU=i' => \$CPU,
             'KEEP_PSLX' => \$KEEP_PSLX,
             "o=s" => \$output_filename,
             );

our $SEE = 0;

$|++;

my $MAX_INTRON = $opt_I || 100000;

my $usage = <<_EOH_;

Script chunks EST alignments into more manageable data sets.

############################# Options ###############################
#
# -g <string>       genomic_seq.db
# -t <string>       transcripts database
# -I <int>          maximum intron length (default: 100000)
# -N <int>          number of top hits (default: $num_top_hits)
# -o <string>       output filename 
#
# --CPU <int>       number of threads (default: 1)
#
# -h this help menu
# -d debug mode
# --KEEP_PSLX      retain the raw blat output files
#
###################### Process Args and Options #####################


_EOH_

    ;


my $genome_db = $opt_g;
my $transcript_db = $opt_t;
my $blat_path = "pblat";
my $util_dir = $FindBin::Bin;

unless ($genome_db && $transcript_db) {
    die "$usage\n";
}

my $output_dir = "pblat_outdir";
if (! -d $output_dir) {
    mkdir($output_dir) or die "Error, cannot mkdir $output_dir";
}


my $checkpt_dir = $output_dir . "/chckpts";
if (! -d $checkpt_dir) {
    mkdir($checkpt_dir) or die "Error, cannot mkdir $checkpt_dir";
}

my $pipeliner = new Pipeliner("-checkpoint_dir" => $checkpt_dir, "-verbose" => 2);


my $ooc_cmd_tmp_out = "tmp-$$-" . int(rand(100000)) . "-out";
my $ooc_cmd = "$blat_path $genome_db $transcript_db -q=rna -dots=100 -maxIntron=$MAX_INTRON -threads=$CPU  -makeOoc=11.ooc $ooc_cmd_tmp_out";
my $ooc_chckpt = "11.ooc.ok";
$pipeliner->add_commands(new Command($ooc_cmd, $ooc_chckpt));

########################
## process pblat search:
########################

my $transcript_file = basename($transcript_db);
my $pslx_file = "$output_dir/$transcript_file.pslx";

my $cmd = "$blat_path $genome_db $transcript_db -q=rna -dots=100 "
    . " -maxIntron=$MAX_INTRON -out=pslx -ooc=11.ooc -threads=$CPU $pslx_file";
        
my $checkpoint_file = basename("$pslx_file.completed");
$pipeliner->add_commands(new Command($cmd, $checkpoint_file));


######################
## get top hits.
######################


my $tophits_file = "$pslx_file.top_${num_top_hits}";
my $tophits_chckpt = basename("$tophits_file.ok");

$cmd = "$util_dir/blat_top_hit_extractor.pl $pslx_file $num_top_hits > $tophits_file";
$pipeliner->add_commands(new Command($cmd, $tophits_chckpt));


$cmd = "$util_dir/pslx_to_gff3.pl < $tophits_file > $tophits_file.gff3";
my $chckpt_file = basename("$tophits_file.gff3.ok");
$pipeliner->add_commands(new Command($cmd, $chckpt_file));


$pipeliner->run();

if($KEEP_PSLX) {
    print STDERR "-retaining $pslx_file\n";
}
else {
    print STDERR "-removing $pslx_file to conserve disk space\n";
    unlink($pslx_file); # these files can be huge. Once have top hits, no longer need all hits (hopefully).
}

$cmd = "ln -sf $tophits_file.gff3 $output_filename";
&Pipeliner::process_cmd($cmd);


print STDERR "done.\n";

exit(0);

