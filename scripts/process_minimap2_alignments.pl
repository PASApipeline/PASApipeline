#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ($FindBin::Bin);
use Pasa_init;

use File::Basename;
use Cwd;

use Pipeliner;

use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);


my $max_intron_length = 100000;

my $usage = <<__EOUSAGE__;

######################################################################
#
#  Required:
#  --genome <string>           target genome to align to
#  --transcripts <string>      cdna sequences to align
#
#  Optional:
#  -I <int>                    maximum intron length (default: $max_intron_length)
#  --gtf <string>              gene structure annotations in gtf format
#  --CPU <int>                 number of threads (default: 2)
#  -o|--output <string>        bam output filename (default: basename(transcripts).mm2.bam)
#
#######################################################################


__EOUSAGE__

    ;

my $genome;
my $transcripts;
my $gtf;
my $CPU = 2;

my $help_flag;
my $output;

&GetOptions( 'h' => \$help_flag,
             'genome=s' => \$genome,
             'transcripts=s' => \$transcripts,
             'gtf=s' => \$gtf,
             'CPU=i' => \$CPU,
             'o|output=s' => \$output,
             'I=i' => \$max_intron_length,
    );


if ($help_flag) {
    die $usage;
}


unless ($genome && $transcripts) {
    die $usage;
}

unless ($output) {
    $output = basename($transcripts) . ".mm2.bam";
}

my $MINIMAP2_CUSTOM_OPTS = $ENV{MINIMAP2_CUSTOM_OPTS} || "";


 main: {
	my $genomeBaseDir = dirname($genome);
	my $genomeName = basename($genome);
	my $genomeDir = "$genomeBaseDir/$genomeName" . ".mm2";
    my $pipeliner = new Pipeliner("-checkpoint_dir" => $genomeDir, "-verbose" => 2);
        
    unless (-d $genomeDir) {
        &process_cmd("mkdir $genomeDir");
    }
    
	my $cwd = cwd();
	
    my $mm2_idx = "$genomeDir/$genomeName.mmi";
	
    my $splice_file = "";
    if ($gtf) {
        $splice_file = "$genomeDir/anno.bed";
    }

    my $minimap2_cmd = "minimap2 -d $mm2_idx -t $CPU $genome";
    my $minimap2_chckpt = "${genomeName}.mmi.ok";

    unless (-e $minimap2_chckpt) {
        $pipeliner->add_commands(new Command($minimap2_cmd, $minimap2_chckpt));
        $pipeliner->run();

        if ($gtf) {
            my $cmd = "paftools.js gff2bed $gtf > $splice_file";
            &process_cmd($cmd);
        }
    }
    
    my $splice_param = "";
    if ($splice_file) {
        $splice_param = "--junc-bed $splice_file";
    }

    my  $genome_samtools_index = "$genomeName.fai";
    
    my $cmd = "bash -c \'set -o pipefail && minimap2 -ax splice $splice_param --secondary=no -O6,24 -B4 -L -t $CPU -cs -ub -G $max_intron_length ${MINIMAP2_CUSTOM_OPTS}  $mm2_idx $transcripts | samtools view -Sbt $genome_samtools_index | samtools sort -o $output && samtools index $output \'";
    &process_cmd($cmd);
    
	exit(0);
}

####
sub process_cmd {
	my ($cmd) = @_;
	
	print STDERR "CMD: $cmd\n";
	#return;

	my $ret = system($cmd);
	if ($ret) {
		die "Error, cmd: $cmd died with ret ($ret)";
	}

	return;
}



