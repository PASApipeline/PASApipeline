#!/usr/bin/env perl

use FindBin;
use lib ("$FindBin::Bin/PerlLib");

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use Gene_obj;
use Nuc_translator;
use Fasta_reader;
use Longest_orf;
use List::Util qw (min max);
use File::Basename;


my $UTIL_DIR = "$FindBin::Bin/util";


my $transcripts_file;

my $min_prot_length = 100;
my $genetic_code;

my $top_ORFs_train = 500;

my $TOP_STRAND_ONLY = 0;

my $help;
my $verbose;
my $search_pfam = "";

my $CPU = 2;
my $RETAIN_LONG_ORFS = 900;


my $usage =  <<_EOH_;

############################# Options #####################################################################
#
# ** Required:
#
# -t <string>                            transcripts.fasta
#
# ** Optional:
#
# -m <int>                               minimum protein length (default: 100)
#
# --search_pfam <string>                 /path/to/pfam_db.hmm to search 
#                                        using hmmscan (which should be accessible via your PATH setting)
#
# -G <string>                            genetic code (default: universal, options: Euplotes, Tetrahymena, Candida, Acetabularia)
#
#
# -h                                     print this option menu and quit
# -v                                     verbose
#
# -S                                     strand-specific (only analyzes top strand)
# -T <int>                               top longest ORFs to train Markov Model (hexamer stats) (default: 500)
#
# --CPU <int>                            number of threads to use (passed on to hmmscan)  (default: 2)
#
# --retain_long_orfs <int>               retain all ORFs found that are of minimum length in nucleotides (default: 900) so 300aa
#
##############################################################################################################

_EOH_

    ;


&GetOptions( 't=s' => \$transcripts_file,
             'm=i' => \$min_prot_length,
             'G=s' => \$genetic_code,
             'h' => \$help,
             'v' => \$verbose,
             'S' => \$TOP_STRAND_ONLY, 
             'T' => \$top_ORFs_train,
             'CPU=i' => \$CPU,
             "search_pfam=s" => \$search_pfam,
             "retain_long_orfs=i" => \$RETAIN_LONG_ORFS,
             );



if ($help) {
    die $usage;
}

if (@ARGV) {
    die "Error, don't understand options: @ARGV";
}



$|++;

our $SEE = $verbose;

unless ($transcripts_file) {
    die "$usage\n";
}

my $hmmscan_prog = "";
if ($search_pfam) {
    
    unless (-s $search_pfam) {
        die "Error, cannot locate pfam database at: $search_pfam";
    }

    ## check for hmmscan utility
    $hmmscan_prog = `which hmmscan`;
    $hmmscan_prog =~ s/\s//g;

    unless ($hmmscan_prog) {
        die "Error, cannot locate 'hmmscan' program in your PATH setting, so can't search pfam";
    }
    
}



my $workdir = "transdecoder.tmp.$$";
mkdir($workdir) or die "Error, cannot mkdir $workdir";

my $NO_REPLACE = 0;

if ($genetic_code) {
    &Nuc_translator::use_specified_genetic_code($genetic_code);
}


my $prefix = "$workdir/longest_orfs";
my $cds_file = "$prefix.cds";
my $gff3_file = "$prefix.gff3";
my $pep_file = "$prefix.pep";


my %orf_lengths;

unless ($NO_REPLACE && -s $pep_file && -s $cds_file && -s $gff3_file) {
	
	open (PEP, ">$pep_file") or die $!;
	open (CDS, ">$cds_file") or die $!; 
	open (GFF, ">$gff3_file") or die $!;
	
	
	my $counter = 0;
	
	my $fasta_reader = new Fasta_reader($transcripts_file);
	while (my $seq_obj = $fasta_reader->next()) {
		
		my $acc = $seq_obj->get_accession();
		my $sequence = $seq_obj->get_sequence();
		
		my $longest_orf_finder = new Longest_orf();
		$longest_orf_finder->allow_5prime_partials();
		$longest_orf_finder->allow_3prime_partials();
		
	    if ($TOP_STRAND_ONLY) {
			$longest_orf_finder->forward_strand_only();
		}
		
		my @orf_structs = $longest_orf_finder->capture_all_ORFs($sequence);
		
		@orf_structs = reverse sort {$a->{length}<=>$b->{length}} @orf_structs;
		
        while (@orf_structs) {
            my $orf = shift @orf_structs;
            
            my $start = $orf->{start};
            my $stop = $orf->{stop};
            
            if ($stop <= 0) { $stop += 3; } # edge issue
            
            my $length = int($orf->{length}/3);
            my $orient = $orf->{orient};
            my $protein = $orf->{protein};
            
            if ($length < $min_prot_length) { next; }
            
            my $coords_href = { $start => $stop };
            
            my $gene_obj = new Gene_obj();
            
            $counter++;
            $gene_obj->populate_gene_object($coords_href, $coords_href);
            $gene_obj->{asmbl_id} = $acc;
            
            my $model_id = "m.$counter";
            my $gene_id = "g.$counter";
            
            
            $gene_obj->{TU_feat_name} = $gene_id;
            $gene_obj->{Model_feat_name} = $model_id;

            
            my $cds = $gene_obj->create_CDS_sequence(\$sequence);
            
            my $got_start = 0;
            my $got_stop = 0;
            if ($protein =~ /^M/) {
                $got_start = 1;
            } 
            if ($protein =~ /\*$/) {
                $got_stop = 1;
            }
            
            my $prot_type = "";
            if ($got_start && $got_stop) {
                $prot_type = "complete";
            } elsif ($got_start) {
                $prot_type = "3prime_partial";
            } elsif ($got_stop) {
                $prot_type = "5prime_partial";
            } else {
                $prot_type = "internal";
            }
            
            $gene_obj->{com_name} = "ORF $gene_id $model_id type:$prot_type len:$length ($orient)";            
            
            print PEP ">$model_id $gene_id type:$prot_type len:$length $acc:$start-$stop($orient)\n$protein\n";
            
            print CDS ">$model_id $gene_id type:$prot_type len:$length\n$cds\n";
            
            print GFF $gene_obj->to_GFF3_format() . "\n";
            

            $orf_lengths{$model_id} = length($cds);
            


        }
	}

    close PEP;
    close CDS;
    close GFF;
    
}

## Train a Markov model based on longest candidate CDS sequences, score all candidates, and select the final set.

# get longest entries
my $top_cds_file = "$cds_file.top_${top_ORFs_train}_longest";
my $cmd = "$UTIL_DIR/get_top_longest_fasta_entries.pl $cds_file $top_ORFs_train > $top_cds_file";
&process_cmd($cmd) unless ($NO_REPLACE && -s $top_cds_file);

#$cmd = "$UTIL_DIR/randomize_sequences.pl $transcripts_file > $transcripts_file.random";
#&process_cmd($cmd) unless ($NO_REPLACE && -s "$transcripts_file.random");

$cmd = "$UTIL_DIR/compute_base_probs.pl $transcripts_file $TOP_STRAND_ONLY > base_freqs.dat";
&process_cmd($cmd);


# get hexamer scores
#$cmd = "$UTIL_DIR/seq_n_background_to_logliklihood_vals.pl $top_cds_file $transcripts_file.random > hexamer.scores";
#&process_cmd($cmd) unless ($NO_REPLACE && -s "hexamer.scores");

$cmd = "$UTIL_DIR/seq_n_baseprobs_to_logliklihood_vals.pl $top_cds_file base_freqs.dat > hexamer.scores";
&process_cmd($cmd);


# score all cds entries
$cmd = "$UTIL_DIR/score_CDS_liklihood_all_6_frames.pl $cds_file hexamer.scores > $cds_file.scores";
&process_cmd($cmd) unless ($NO_REPLACE && -s "$cds_file.scores");

# run pfam
my %has_pfam_hit;
if ($search_pfam) {
    my $cmd = "$hmmscan_prog --cpu $CPU --noali --cut_nc --acc --notextw --tblout $pep_file.pfam.dat $search_pfam $pep_file ";
    
    &process_cmd($cmd);

    open (my $fh, "$pep_file.pfam.dat") or die "Error, cannot open file: $pep_file.pfam.dat";
    while (<$fh>) {
        chomp;
        my @x = split(/\s+/);
        my $orf_acc = $x[2];
        $has_pfam_hit{$orf_acc} = 1;
    }
    close $fh;
}


# get accs for best entries
my $acc_file = "$cds_file.scores.selected";
{
	open (my $ofh, ">$acc_file") or die "Error, cannot write to $acc_file";
	open (my $ifh, "$cds_file.scores") or die "Error, cannot open file $cds_file.scores";
	while (<$ifh>) {
		chomp;
		my ($acc, @scores) = split(/\t/);
		
		my $score_1 = shift @scores;
		my $max_score_other_frame = max(@scores);
		if ($has_pfam_hit{$acc} 
            || 
            $orf_lengths{$acc} >= $RETAIN_LONG_ORFS
            ||
            ($score_1 > 0 && $score_1 > $max_score_other_frame)
            ) { 
			print $ofh "$acc\n";
		}
	}
	close $ifh;
	close $ofh;
}

# index the current gff file:
$cmd = "$UTIL_DIR/index_gff3_files_by_isoform.pl $gff3_file";
&process_cmd($cmd);

# retrieve the best entries:
$cmd = "$UTIL_DIR/gene_list_to_gff.pl $acc_file $gff3_file.inx > $cds_file.best_candidates.gff3";
&process_cmd($cmd);

{
    my $final_output_prefix = basename($transcripts_file) . ".transdecoder";
    
    # exclude shadow orfs (smaller orfs in different reading frame that are eclipsed by longer orfs)
    $cmd = "$UTIL_DIR/remove_eclipsed_ORFs.pl $cds_file.best_candidates.gff3 > $final_output_prefix.gff3";
    &process_cmd($cmd);
    


    ## write final outputs:
    
    ## make a BED file for viewing in IGV
    my $gff3_file = "$final_output_prefix.gff3";
    my $bed_file = $gff3_file;
    $bed_file =~ s/\.gff3$/\.bed/;
    $cmd = "$UTIL_DIR/gff3_file_to_bed.pl $gff3_file > $bed_file";
    &process_cmd($cmd);
    
    
    # make a peptide file:
    my $best_pep_file = $gff3_file;
    $best_pep_file =~ s/\.gff3$/\.pep/;
    $cmd = "$UTIL_DIR/gff3_file_to_proteins.pl $gff3_file $transcripts_file > $best_pep_file";
    &process_cmd($cmd);



    # make a CDS file:
    my $best_cds_file = $best_pep_file;
    $best_cds_file =~ s/\.pep$/\.cds/;
    $cmd = "$UTIL_DIR/gff3_file_to_proteins.pl $gff3_file $transcripts_file CDS > $best_cds_file";
    &process_cmd($cmd);
    
}

print STDERR "transdecoder is finished.\n";


exit(0);


####
sub process_cmd {
	my ($cmd) = @_;

	print "CMD: $cmd\n";
	my $ret = system($cmd);

	if ($ret) {
		die "Error, cmd: $cmd died with ret $ret";
	}
	
	return;

}





