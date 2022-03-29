#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use ConfigFileReader;
use File::Basename;
use Pipeliner;

use Getopt::Long qw(:config no_ignore_case bundling pass_through);

my $usage = <<__EOUSAGE__;

####################################################################################
#
#
#  --align_assembly_config <string>  align_assembly.config file
#
#  --annot_compare_config <string>   annot_compare.config file
#
# ## all optional:
#
# --ALIGNERS <string>   "blat", "gmap", or "minimap2" (default: blat,gmap,minimap2)
#
# --CPU <int>      number of threads to use (default: 2)
#
# --TRANSDECODER   run TRANSDECODER to identify candidate full-length transcripts
#
# --just_align_assembly    only run the initial alignment assembly portion of the pipeline.
#
# -s <int>                 resume alignAssembly run at step
#
# -N <int>              number of top scoring spliced alignments (default: 1)
#
# --stringent_alignment_overlap <int>     (suggested: 30.0)  overlapping transcripts must have this min % overlap to be clustered.
# --gene_overlap <int>     (suggested: 50.0)  transcripts overlapping existing gene annotations are clustered.  Intergenic alignments are clustered by default mechanism.                
#
####################################################################################


__EOUSAGE__

    ;


my $help_flag;

my $ALIGNERS = "blat,gmap,minimap2";
my $CPU = 4;
my $TRANSDECODER;
my $JUST_ALIGN_ASSEMBLY;
my $resume_step;
my $num_top_hits = 1;
my $stringent_alignment_overlap;
my $gene_overlap;
my $align_assembly_config_file;
my $annot_compare_config_file;


&GetOptions ( 'h' => \$help_flag,
              'align_assembly_config=s' => \$align_assembly_config_file,
              'annot_compare_config=s' => \$annot_compare_config_file,
                            
              'CPU=i' => \$CPU,
              'TRANSDECODER' => \$TRANSDECODER,
              'ALIGNERS=s' => \$ALIGNERS,
              'just_align_assembly' => \$JUST_ALIGN_ASSEMBLY,
              's=i' => \$resume_step,
              'N=i' => \$num_top_hits,
              'stringent_alignment_overlap=i' => \$stringent_alignment_overlap,
              'gene_overlap=i' => \$gene_overlap,
              

    );


if ($help_flag) {
    die $usage;
}

unless ($align_assembly_config_file && $annot_compare_config_file) {
    die $usage;
}


my %config = &readConfig($align_assembly_config_file);
my $DBname = $config{DATABASE} or die "Error, couldn't extract DATABASE from config file: $align_assembly_config_file";

$DBname = basename($DBname);

 main: {

     my $checkpoints_dir = "__chkpts_$DBname";
     my $pipeliner = new Pipeliner(-verbose=>2,
                                   -checkpoint_dir=>$checkpoints_dir,
                                   -cmds_log=> "$checkpoints_dir.cmds_log",
         );
     
     
     if (-s "genome_sample.fasta.gz" && ! -s "genome_sample.fasta") {
         &process_cmd("gunzip -c genome_sample.fasta.gz > genome_sample.fasta");
     }
     
    
  # goto annot_compare;
  #  goto annot_compare_R2;
  
   alignment_assembly: 
     {
		
		print "********* Running Alignment Assembly ************\n";
		
		# "-C -r" will drop db if exists
		my $cmd = "../Launch_PASA_pipeline.pl -c $align_assembly_config_file -C -r -R -g genome_sample.fasta -t all_transcripts.fasta.clean -T -u all_transcripts.fasta -f FL_accs.txt --ALIGNERS $ALIGNERS --CPU $CPU -N $num_top_hits --TDN tdn.accs  --IMPORT_CUSTOM_ALIGNMENTS_GFF3 custom_alignments.gff3 ";
		
        if ($TRANSDECODER) {
            $cmd .= " --TRANSDECODER ";
        }

        if ($resume_step) {
            $cmd .= " -s $resume_step ";
        }

        if ($stringent_alignment_overlap) {
            $cmd .= " --stringent_alignment_overlap $stringent_alignment_overlap ";
        }
        elsif ($gene_overlap) {
            $cmd .= " --gene_overlap $gene_overlap --annots orig_annotations_sample.gff3 ";
        }
        
		$pipeliner->add_commands(new Command($cmd, "align_assembly.ok"));
	}

     $pipeliner->run();
     
     
    
  comprehensive_transcriptome_build:
     {
         print "********** Building comprehensive transcriptome database ***********\n";
         my $cmd = "../scripts/build_comprehensive_transcriptome.dbi -c $align_assembly_config_file -t all_transcripts.fasta.clean";
         $pipeliner->add_commands(new Command($cmd, "build_compreh.ok"));
     }

     $pipeliner->run();
     
     if ($JUST_ALIGN_ASSEMBLY) {
         print STDERR "-stopping now, after alignment assembly.\n";
         exit(0);
     }
     


	
   annot_compare_R1:
     { ## Annotation comparisons:
         
         ## First round, using pre-existing gene structure annotations:
         
         print "******** Comparing Annotations to Alignment Assemblies ***********\n";
         my $cmd = "../Launch_PASA_pipeline.pl -c $annot_compare_config_file -g genome_sample.fasta -t all_transcripts.fasta.clean -A -L --annots orig_annotations_sample.gff3 --CPU $CPU";
         $pipeliner->add_commands(new Command($cmd, "annot_compare_R1.ok"));
         
     }
	
     $pipeliner->run();
     
   annot_compare_R2:
     {
         
         ## Run it again, using the output of updated genes from the first round.
         ## maybe capture a few more updates, and at the very least, verify that the inital update worked!
         
         
         print "********** Loading Updated Gene Annotations ************\n";
         # get the file containing the updates:
         my @update_files = grep { /gene_structures_post_PASA_updates.\d+.gff3/ } `ls -t`;
         chomp @update_files;
         my $recent_update_file = shift @update_files;
         unless ($recent_update_file) { 
             die "Error, couldn't identify the gff3 file containing the pasa-based annotation updates!"; 
         }
         
         print "******** Comparing Annotations to Alignment Assemblies ***********\n";
         my $cmd = "../Launch_PASA_pipeline.pl -c $annot_compare_config_file -g genome_sample.fasta -t all_transcripts.fasta.clean -A -L --annots $recent_update_file --CPU $CPU ";

         $pipeliner->add_commands(new Command($cmd, "annot_compare_R2.ok"));
     }
    
     $pipeliner->run();
     
   alt_splice_analysis:
     
     {
         print "*********** Running Analysis of Alternative Splicing *******\n";
         my $cmd = "../Launch_PASA_pipeline.pl -c $annot_compare_config_file -g genome_sample.fasta -t all_transcripts.fasta.clean --CPU $CPU --ALT_SPLICE";
         $pipeliner->add_commands(new Command($cmd, "alt_splicing.ok"));
         
     }

     $pipeliner->run();
     
     
   find_orfs_in_pasa_assemblies:
     
     {
         
         print "***********  Finding ORFs in PASA assemblies **************\n";
         my $cmd = "../scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta $DBname.assemblies.fasta --pasa_transcripts_gff3 $DBname.pasa_assemblies.gff3";
         $pipeliner->add_commands(new Command($cmd, "find_orfs_in_pasa.ok"));
         
     }

     $pipeliner->run();
     
     
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
