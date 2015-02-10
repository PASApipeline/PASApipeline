#!/usr/bin/env perl

use strict;
use warnings;
use Carp;


## Run script with no parameters to demonstrate the default pipeline.
## To include the sim4 chaser (slow, capture a few extra valid alignments), include any parameter value.

my $run_sim4_chaser = $ARGV[0];


main: {
 
  # goto annot_compare;
  
  initial_cleanup:
	{
		## Purge the current sample mysql database if it exists from a previous run of this pipeline. Start fresh.
		
		my $cmd = "../scripts/drop_mysql_db_if_exists.dbi -c alignAssembly.config";
		&process_cmd($cmd);
	}
	
	
  alignment_assembly: 
	{
		
		print "********* Running Alignment Assembly ************\n";
		
		
		my $cmd = "../scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g genome_sample.fasta -t all_transcripts.fasta.clean -T -u all_transcripts.fasta -f FL_accs.txt ";
		
		if (defined $run_sim4_chaser) {
			$cmd .= " --APPLY_SIM4_CHASER ";
		}
		
		&process_cmd($cmd);
	}
	
  annot_compare_R1:
	{ ## Annotation comparisons:
		
		## First round, using pre-existing gene structure annotations:
		
		print "******** Comparing Annotations to Alignment Assemblies ***********\n";
		my $cmd = "../scripts/Launch_PASA_pipeline.pl -c annotCompare.config -g genome_sample.fasta -t all_transcripts.fasta.clean -A -L --annots_gff3 orig_annotations_sample.gff3 ";
		&process_cmd($cmd);
		
	}
	
	
	
  annot_compare_R2:
	{
		
		## Run it again, using the output of updated genes from the first round.
		## maybe capture a few more updates, and at the very least, verify that the inital update worked!
		
		
		print "********** Loading Updated Gene Annotations ************\n";
		# get the file containing the updates:
		my @update_files = grep { /gene_structures_post_PASA_updates.\d+.gff3/ } `ls -t`;
		my $recent_update_file = shift @update_files;
		unless ($recent_update_file) { 
			die "Error, couldn't identify the gff3 file containing the pasa-based annotation updates!"; 
		}
		
		print "******** Comparing Annotations to Alignment Assemblies ***********\n";
		my $cmd = "../scripts/Launch_PASA_pipeline.pl -c annotCompare.config -g genome_sample.fasta -t all_transcripts.fasta.clean -A -L --annots_gff3 $recent_update_file";
		&process_cmd($cmd);
	}
		
	
	exit(0);
}



####
sub process_cmd {
    my ($cmd) = @_;
    
    print "CMD: $cmd\n\n";
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd: $cmd died with ret($ret) ";
    }

    return;
}
