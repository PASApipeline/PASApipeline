#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use Carp;
use Data::Dumper;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
#use lib ("$ENV{EUK_MODULES}");
use Fasta_reader;
use GFF3_alignment_utils;

my $usage = <<__EOUSAGE__;

##############################################################################################
#
#  Details on input transcripts:
#
#  --trans_fasta <string>              input transcripts in fasta format.
#  --trans_validations <string>        gmap validations file 
#  --trans_blastx_m8 <string>          swissprot blastx m8 file
#  --trans_TD_pep <string>             transdecoder pep file for input transcripts
#  --trans_chim_jaccard_summary <string>        jaccard stats for chims
#  --trans_cov_info <string>           read coverage info
#
#  
# TDN details:
#  --RSEM_for_TDN <string>             RSEM output file for TDN only

# PASA details
#
#  --pasa_gff3 <string>                pasa_db.pasa_assemblies.denovo_transcript_isoforms.gff3 file
#  --pasa_TD_pep <string>              pasa transdecoder pep file (including the '>PASA|' prefix)
#
# Genome details:
#  --genome_gaps <string>              positions of sequencing gaps in the genome.
#  --scaffold_lengths <string>         lengths of scaffolds file. Format:  scaffold(tab)length
#
#  Details for comparisons between PASA assemblies and Input transcripts
#
#  --cdhit_clstr <string>              cdhit results
#                                         1. make pasa.transdecoder.pep have accessions with '>PASA|' prefix.
#                                         2. combine pasa.transdecoder.pep with input_trans.transdecoder.pep > both.PEP
#                                         cdhit run like so: "cd-hit -i both.PEP -c 0.98 -o cdhit -p 1"
#
#
#  Optional:
#
#  --proteomics_files <string>         comma delimited list of files containing proteomic hit information
#
#  -F                                  report failed alignments only
#
##################################################################################################


__EOUSAGE__

    ;


my $help_flag;

my $trans_fasta_file;
my $trans_validations_file;
my $trans_blastx_m8_file;
my $trans_TD_pep_file;
my $trans_chim_jaccard_summary;
my $trans_cov_info_file;

my $RSEM_TDN_file;

my $pasa_gff3_file;
my $pasa_TD_pep_file;

my $genome_gaps_file;
my $scaffold_lengths_file;

my $cdhit_clstr_file;

my $REPORT_FAILURES_ONLY_FLAG = 0;

my $proteomics_files = "";

## Gap info
my $MIN_GAP_LEN = 50;
my $DIST_TO_GAP_REPORT = 100;
my $DIST_TO_SCAFF_EDGE_REPORT = 1000;


## Good chimera candidates require:
my $MIN_PAIR_SPAN = 10;
my $MIN_JACCARD = 0.05;
my $MIN_EXPR_QUANTILE = 0.1;
my $MIN_ALIGN_PER_ID = 98;

## Proteomics
my $MIN_TAGS_PROTEOMICS_REPORT = 2;



&GetOptions ( 'h' => \$help_flag,
              
              'trans_fasta=s' => \$trans_fasta_file,
              'trans_validations=s' => \$trans_validations_file,
              'trans_blastx_m8=s' => \$trans_blastx_m8_file,
              'trans_TD_pep=s' => \$trans_TD_pep_file,
              'trans_chim_jaccard_summary=s' => \$trans_chim_jaccard_summary,
              'trans_cov_info=s' => \$trans_cov_info_file,
              
              'RSEM_for_TDN=s' => \$RSEM_TDN_file,

              'genome_gaps=s' => \$genome_gaps_file,
              'scaffold_lengths=s' => \$scaffold_lengths_file,

              'pasa_gff3=s' => \$pasa_gff3_file,
              'pasa_TD_pep=s' => \$pasa_TD_pep_file,
              
              'cdhit_clstr=s' => \$cdhit_clstr_file,
              
              'F' => \$REPORT_FAILURES_ONLY_FLAG,
              
              'proteomics_files=s' => \$proteomics_files,
              
              );


if ($help_flag) {
    die $usage;
}

unless ($trans_fasta_file
        && $trans_validations_file 
        && $trans_blastx_m8_file
        && $trans_TD_pep_file
        && $trans_chim_jaccard_summary
        && $trans_cov_info_file

        # && $RSEM_TDN_file

        && $genome_gaps_file
        && $scaffold_lengths_file

        && $pasa_gff3_file
        && $pasa_TD_pep_file
        
        && $cdhit_clstr_file
        ) {

    die $usage;
}


## ensure we can find all input files specified.
foreach my $file ($trans_fasta_file, $trans_validations_file, $trans_blastx_m8_file,
                  $trans_TD_pep_file, $trans_chim_jaccard_summary, $trans_cov_info_file,
                  
                  $RSEM_TDN_file,

                  $pasa_gff3_file, $pasa_TD_pep_file,
                  
                  $genome_gaps_file,

                  $cdhit_clstr_file) {

    unless (-s $file) {
        die "Error, cannot find file $file";
    }
}


main: {

    my %transcript_to_info;
    
    print STDERR "-populating transcript info from fasta file: $trans_fasta_file\n";
    &populate_transcript_base_info(\%transcript_to_info, $trans_fasta_file);

    print STDERR "-populating pasa assembly info from gff3 file: $pasa_gff3_file\n";
    &populate_pasa_structure_info(\%transcript_to_info, $pasa_gff3_file);
    
    if ($proteomics_files) {
        &populate_proteomics_info(\%transcript_to_info, $proteomics_files);
    }
 

    print STDERR "-populating_chimera_jaccard_info\n";
    &populate_chimera_jaccard_info(\%transcript_to_info, $trans_chim_jaccard_summary);
    
    
    my %transdecoder_acc_to_transcript_acc;
    print STDERR "-populating transdecoder results for: $trans_TD_pep_file\n";
    &populate_transdecoder_orfs(\%transcript_to_info, \%transdecoder_acc_to_transcript_acc, $trans_TD_pep_file);

    print STDERR "-populating transdecoder results for: $pasa_TD_pep_file\n";
    &populate_transdecoder_orfs(\%transcript_to_info, \%transdecoder_acc_to_transcript_acc, $pasa_TD_pep_file);

    
    print STDERR "-populating alignment and validation status info from $trans_validations_file\n";
    &populate_alignment_and_validation_info(\%transcript_to_info, $trans_validations_file);
    
    print STDERR "-populating read coverage info from $trans_cov_info_file\n";
    &populate_read_coverage_info(\%transcript_to_info, $trans_cov_info_file);

    if ($RSEM_TDN_file) {
        print STDERR "-populating RSEM for TDN transcripts\n";
        &populate_RSEM_for_TDN(\%transcript_to_info, $RSEM_TDN_file);
    }
    
    print STDERR "-reading homology info from blastx file: $trans_blastx_m8_file";
    &populate_blastx_info(\%transcript_to_info, $trans_blastx_m8_file);

    print STDERR "-reading scaffold lengths file\n";
    my %scaffold_lengths = &populate_scaffold_lengths($scaffold_lengths_file);
    

    print STDERR "-reading genome sequencing gaps file\n";
    my %scaffold_to_gap_positions = &populate_genome_seq_gaps($genome_gaps_file);

    
    my %cdhit_clusters;
    print STDERR "-parsing cdhit protein clustering results: $cdhit_clstr_file\n";
    &populate_cdhit_clustering(\%transcript_to_info, \%transdecoder_acc_to_transcript_acc, \%cdhit_clusters, $cdhit_clstr_file);


    print STDERR "-now analyzing data.\n\n";
    
    # print header:
    print join("\t", "#class_token", "trans_acc", "median_read_cov", "quintile", "num_aligns", "blastx_hit", "has_valid_alignment", "alignments...") . "\n";
    

    ## summarize transcripts.
    foreach my $struct (values %transcript_to_info) {
        
        unless ($struct->{is_transcript}) { next; }
        
        my $acc = $struct->{acc};
        my $median_cov = $struct->{median_read_cov};
        my $quintile_read_cov = sprintf("%.3f", $struct->{quintile_read_cov});

        my $blastx_hit_info = $struct->{blastx_hit} || ".";
        
        my $mapped = $struct->{mapped};
        my $valid = ($struct->{valid}) ? "valid:YES": "valid:NO";
        
        if ($REPORT_FAILURES_ONLY_FLAG && $valid eq "valid:YES") { next; }
        
        my @summary_info = (".", $acc, $median_cov, $quintile_read_cov, $mapped, $blastx_hit_info, $valid);
        
        if ($valid eq "valid:YES") {
            $summary_info[0] = "OK";
        }
        
        if ($mapped) { 
            
            my @alignments = @{$struct->{alignment_structs}};
        
            my @alignment_texts;
            my @alignment_per_ids; # used in the chimera check below.
            foreach my $alignment_struct (@alignments) {
            
                my $mapping_info = $alignment_struct->{scaffold} . ":" 
                    . $alignment_struct->{lend} . "-" . $alignment_struct->{rend}
                    . " " . $alignment_struct->{alignment}
                . " PerID: " . $alignment_struct->{per_id}
                . " PerAlign: " . $alignment_struct->{per_aligned};
                
                push (@alignment_per_ids, $alignment_struct->{per_id});
                
                if (my $failed_alignment_comment = $alignment_struct->{failed_alignment_comment}) {
                    $mapping_info .= " ** $failed_alignment_comment ** ";
                }
                                
                push (@alignment_texts, $mapping_info);
                
            }

            my $alignment_text = join(" |||| ", @alignment_texts);
            
            push (@summary_info, $alignment_text);
            

            if (1) { #$valid eq "valid:NO") {
                ## Options to consider:
                

                ## all based on the transdecoder orf clustering

                my %pasa_cdhit_compare = &compare_to_PASA_cdhit_orfs($struct, \%transcript_to_info, \%transdecoder_acc_to_transcript_acc, \%cdhit_clusters);

                ## check for chimera
                if ($mapped > 1) {

                    $summary_info[0] .= "|CHIMERIC_ALIGNMENT";
                    
                    my $chimera_info_href = $struct->{chimera_info};
                    
                    my $jaccard_val = $chimera_info_href->{jaccard};
                    my $both_count = $chimera_info_href->{both_count};
                    
                    push (@summary_info, "ChimeraInfo: JaccardVal=$jaccard_val, PairSpanCount=$both_count");
                    
                    ## check to see if it's a good chimera candidate.
                    if ($jaccard_val >= $MIN_JACCARD
                        && $both_count >= $MIN_PAIR_SPAN
                        && $quintile_read_cov >= $MIN_EXPR_QUANTILE
                        && (! &has_overlapping_alignments($struct))
                        
                        ) {
                        
                        unless (grep { $_ < $MIN_ALIGN_PER_ID } @alignment_per_ids) {
                            
                            $summary_info[0] .= "|CHIMERA_Candidate";
                            
                        }
                    }
                    
                    
                    if (my $chimera_info = &chimeric_transcript_breaks_orf($struct, \%transcript_to_info, \%transdecoder_acc_to_transcript_acc) ) {
                    
                        ## Include multiple categories here
                        
                        # chimeric transcript, likely a transcript misassembly
                        # chimeric transcript, well supported by read pairings
                        #     -chimera split orf
                        #     -doesn't interrupt orf
                        
                        $summary_info[0] .= "|CHIMERA_SPLIT_ORF";
                        push (@summary_info, $chimera_info);
                    
                    }
                    
                    

                }

                                
                if (my $near_gap_note = &near_sequencing_gap($struct, \%scaffold_to_gap_positions)) {
                    
                    $summary_info[0] .= "|Adj_to_GAP";
                    push (@summary_info, $near_gap_note);
                    
                }
                
                if (my @near_edge_notes = &near_scaffold_edge($struct, \%scaffold_lengths)) {
                    
                    #print STDERR Dumper(\@near_edge_notes);
                    
                    $summary_info[0] .= "|Adj_to_Scaffold_Edge";
                    push (@summary_info, @near_edge_notes);
                    
                    my %scaffolds_aligned;
                    foreach my $note (@near_edge_notes) {
                        if ($note =~ /^(alignment:[^:]+)/) {
                            my $scaff = $1;
                            $scaffolds_aligned{$scaff}++;
                        }
                        else {
                            die "Error, cannot parse scaffold info from note: $note";
                        }
                    }
                    

                    if (scalar(keys %scaffolds_aligned) > 1) {
                        $summary_info[0] .= "|Scaffold_Bridge";
                    }
                    
                }
                
                
                if (my $note = &ORF_partially_aligned($struct, \%transcript_to_info, \%transdecoder_acc_to_transcript_acc)) {
                    
                    $summary_info[0] .= "|ORF_partially_aligned";
                    push (@summary_info, $note);
                }
                
                ## Examine the type of alignment error
                my @alignments = @{$struct->{alignment_structs}};
                
                my @tokens;
                
                foreach my $alignment (@alignments) {
                    
                    my $failure_comment = $alignment->{failed_alignment_comment} || "";
                    if ($failure_comment  =~ /Only .* identity on average/) {
                        push (@tokens, "|Weak_alignment");
                    }
                    elsif ($failure_comment =~ /Only .* is aligned, less/) {
                        push (@tokens, "|Partial_alignment");
                    }
                    elsif ($failure_comment =~ /Splice site validations failed|splice boundary position/) {
                        push (@tokens, "|Splice_problem");
                    }
                    elsif ($failure_comment =~ /Intron/) {
                        push (@tokens, "|Long_intron");
                    }
                    elsif ($failure_comment =~ /Incontiguous alignment/) {
                        push (@tokens, "|Incontiguous_alignment");
                    }
                }
                
                $summary_info[0] .= join("", @tokens);
                

                ## Compare to any existing cdhit-clustered PASA ORF
                                
                ## same coding content as the PASA assembly, so presumably fewer worries
                if ($pasa_cdhit_compare{same_as_pasa}) {
                    
                    $summary_info[0] .= "|PASA_SAME";
                    push (@summary_info, "no better than the current PASA ORF");
                }
                                
                # - mapped location already represented by a valid PASA assembly (is de novo better than pasa?)
                #       -encapsulates a cdhit pasa orf from a pasa assembly mapped at the same genomic location.
                elsif (my $per_length = $pasa_cdhit_compare{longer_than_pasa}) { 
                    
                    my $note = $pasa_cdhit_compare{longer_than_pasa_note};

                    $summary_info[0] .= "|TRANS_better_than_PASA";
                    push (@summary_info, "trans-encoded orf is $per_length\% longer than the PASA orf, $note");
                }
                
                ## a PASA assembly ORF already represents this ORF.

                elsif ($per_length = $pasa_cdhit_compare{shorter_than_pasa}) { 
                    
                    my $note = $pasa_cdhit_compare{shorter_than_pasa_note};

                    $summary_info[0] .= "|TRANS_inferior_to_PASA";
                    push (@summary_info, "trans-encoded orf is only $per_length\% of the PASA orf length, $note");
                }
                

                # - mapped location not represented by a PASA assembly-containing ORF (does genome interfere w/ transcript representation?)
                #       
                #       does it have an ORF?
                #
                #       
                
                if (%{$struct->{transdecoder_orfs}}) {
                    
                    $summary_info[0] .= "|Mapped_w_ORF";
                    
                    my $transdecoder_orf_summary = &summarize_orfs_found(values %{$struct->{transdecoder_orfs}});
                    
                    push (@summary_info, $transdecoder_orf_summary);
                    
                }
                else {
                    
                    $summary_info[0] .= "|Mapped_w/o_ORF";
                }
                                
            }

        }
        else {

            ##############################
            ## not mapped to the genome ##
            ##############################

            # Is it a 'missing' gene?
            
            if (%{$struct->{transdecoder_orfs}}) {
                
                $summary_info[0] .= "|MISSING_w_ORF";
                
                my $transdecoder_orf_summary = &summarize_orfs_found(values %{$struct->{transdecoder_orfs}});
                
                push (@summary_info, $transdecoder_orf_summary);
                
            }
            else {
                
                $summary_info[0] .= "|MISSING_w/o_ORF";
            }
            
        }

        ## Check proteomics data.
        if (my %proteomics_info = %{$struct->{proteomics_info}}) {
            
            my $got_proteomics_hit_plus = 0;
            my $got_proteomics_hit_minus = 0;
            
            foreach my $sample_type (keys %proteomics_info) {
                
                my $hits = $proteomics_info{$sample_type};
                push (@summary_info, "proteomics($sample_type)=$hits");
                
                my ($plus_hits, $minus_hits) = split(/;/, $hits);
                
                my @p_hits = split(/,/, $plus_hits);
                my @m_hits = split(/,/, $minus_hits);
                
                if (scalar @p_hits >= $MIN_TAGS_PROTEOMICS_REPORT) {
                    $got_proteomics_hit_plus = 1;
                }
                if (scalar @m_hits >= $MIN_TAGS_PROTEOMICS_REPORT) {
                    $got_proteomics_hit_minus = 1;
                }
            }
            
            if ($got_proteomics_hit_plus) {
                $summary_info[0] .= "|PROTEO_PLUS";
            }
            if ($got_proteomics_hit_minus) {
                $summary_info[0] .= "|PROTEO_MINUS";
            }
        }
        

        ### Finally, report result
        
        print join("\t", @summary_info) . "\n";
        
    }
    
    
    

    exit(0);
}

####
sub populate_transcript_base_info {
    my ($transcript_to_info_href, $trans_fasta_file) = @_;

    my $fasta_reader = new Fasta_reader($trans_fasta_file);
    
    while (my $seq_obj = $fasta_reader->next()) {
        
        my $acc = $seq_obj->get_accession();
        my $sequence = $seq_obj->get_sequence();

        $transcript_to_info_href->{$acc} = { 

            is_transcript => 1, # pasa assemblies are set to zero
            
            acc => $acc,
            
            length => length($sequence),
            

            # genome alignment info
            mapped => 0,  # store number of alignments
            valid => 0,
            alignment_structs => [], # has format {alignment => undef,
            #                                          valid => 0,
            #                                       scaffold => undef,
            #                                           lend => undef,
            #                                           rend => undef,
            #                                          mlend => undef,
            #                                          mrend => undef,
            #                                   },
            
            # read coverage
            median_read_cov => 0,
            quintile_read_cov => 0,

            # blastx homology
            blastx_hit => "",

            cdhit_encapsulates_shorter => [],  # list of orf structs
                                               # { acc => str, 
            #                                      aa_len => int }
            cdhit_encapsulated_by_longer => [],
            cdhit_same_length => [],
            
            transdecoder_orfs => {},  # key on acc, # orf structs
            #                               { acc => str,
            #                                 aa_len => int,
            #                                 lend   => int,  # coordinate in transcript space
            #                                 rend   => int,  
            #                                 orient => [+-],
            #                                 cdhit_clust_id => string
            #                               }
            

            chimera_info => {},


            proteomics_info => {},
            
        };

    }

    return;
}


####
sub populate_alignment_and_validation_info {
    my ($transcript_to_info_href, $trans_validations_file) = @_;
    
    
    open (my $fh, $trans_validations_file) or die "Error, cannot open file $trans_validations_file";
    while (<$fh>) {
        unless (/\w/) { next; }
        if (/^\#/) { next; }
        chomp;
        my @x = split(/\t/);
        my $acc = $x[1];
        
        my $scaffold = $x[4];
        my $coord_span = $x[9];
        my ($lend, $rend) = split(/-/, $coord_span);

        my $avg_per_id = $x[10];
        my $per_aligned = $x[11];
        
        my $valid_status = $x[8];
        
        my $alignment = $x[12];
        my ($mlend, $mrend) = &parse_cdna_span_from_alignment_text($alignment);


        my $comment = $x[13];
        
        $acc =~ s/\.path\d+$//; ## accommodate multiple paths here
        
        my $struct = $transcript_to_info_href->{$acc} or die "Error, no struct for acc: $acc ";

        $struct->{mapped}++;
        if ($valid_status eq "OK") {
            $struct->{valid} = 1;
        }
        my $alignment_struct = { alignment => $alignment,
                                 valid => ($valid_status eq "OK") ? 1:0,
                                 scaffold => $scaffold,
                                 lend => $lend,
                                 rend => $rend,
                                 mlend => $mlend,
                                 mrend => $mrend,
                                 failed_alignment_comment => $comment,
                             
                                 per_id => sprintf("%.2f", $avg_per_id),
                                 per_aligned => sprintf("%.2f", $per_aligned),
                             };
        
        
        push (@{$struct->{alignment_structs}}, $alignment_struct);
        
    }

    return;
}


####
sub populate_read_coverage_info {
    my ($transcript_to_info_href, $trans_cov_info_file) = @_;

    open (my $fh, $trans_cov_info_file) or die $!;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $acc = $x[0];
        my $median_read_cov = $x[2];
        
        my $struct = $transcript_to_info_href->{$acc} or confess "Error, no struct for acc: $acc ";
        
        $struct->{median_read_cov} = $median_read_cov;

    }
    close $fh;

    &set_expression_quantiles($transcript_to_info_href);

    return;
}

####
sub populate_RSEM_for_TDN {
    my ($transcript_to_info_href, $RSEM_TDN_file) = @_;

    open (my $fh, $RSEM_TDN_file) or die $!;
    my $header = <$fh>;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $acc = $x[0];
        unless ($acc =~ /^DN\|/) { next; }
        my $count = $x[4];
        my $length = $x[2];

        my $mean_count = sprintf("%.3f", ($count*25)/($length-24)); # kmer coverage
        
        my $struct = $transcript_to_info_href->{$acc} or confess "Error, no struct for acc: $acc ";
        
        $struct->{median_read_cov} = $mean_count; # not the median, but reusing the variable from earlier.
        
    }
    close $fh;

    &set_expression_quantiles($transcript_to_info_href);
    
    return;
}

####
sub set_expression_quantiles {
    my ($transcript_to_info_href) = @_;
    
    ## determine quintiles
    my @structs = values %$transcript_to_info_href;
    @structs = grep { $_->{is_transcript} == 1 } @structs;  ## dangerous grouping together the pasa assemblies and alignments this way.
    
    ## initialize all
    foreach my $struct (@structs) {
        $struct->{quintile_read_cov} = 0;
    }
    

    @structs = sort {$a->{median_read_cov}<=>$b->{median_read_cov}} grep {$_->{median_read_cov} >= 1 } @structs;
    
    my $num_structs = scalar(@structs);

    for (my $i = 0; $i <= $#structs; $i++) {
        
        $structs[$i]->{quintile_read_cov} = $i/$num_structs;
    }
    

    return;
}
    
####
sub populate_blastx_info {
    my ($transcript_to_info_href, $trans_blastx_m8_file) = @_;

    open (my $fh, $trans_blastx_m8_file) or die $!;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $acc = $x[0];
        my $hit = $x[1];
        my $per_id = $x[2];
        my $evalue = $x[10];
        my $descr = $x[12];

        if ($descr) {
            $hit .= " ($descr)";
        }
        
        my $struct = $transcript_to_info_href->{$acc} or confess "Error, no struct for acc: $acc";
        
        $struct->{blastx_hit} = join("||", $hit, $per_id, $evalue);

    }

    close $fh;

    return;
}

####
sub populate_cdhit_clustering {
    my ($transcript_to_info_href, $transdecoder_acc_to_transcript_acc_href, $cdhit_clusters_href, $cdhit_clstr_file) = @_;

    my @entries;
    

    my $cluster_counter = 0;
    open (my $fh, $cdhit_clstr_file) or die $!;
    while (<$fh>) {
        chomp;
        if (/^>Cluster/) {
            if (@entries) {
                $cluster_counter++;
                &_process_cdhit_entries($transcript_to_info_href, $transdecoder_acc_to_transcript_acc_href, \@entries, $cdhit_clusters_href, $cluster_counter);
            }
            @entries = ();
        }
        else {
            my @x = split(/\s+/);
            my $len = $x[1];
            $len =~ s/aa,$//;
            my $acc = $x[2];
            $acc =~ s/>//;
            my $align_descr = pop @x;
            
            my $entry = { aa_len => $len,
                          acc => $acc,
                      };
            push (@entries, $entry);
            
        }
    }
    close $fh;

    # get last ones
    if (@entries) {
        $cluster_counter++;
        &_process_cdhit_entries($transcript_to_info_href, $transdecoder_acc_to_transcript_acc_href, \@entries, $cdhit_clusters_href, $cluster_counter);
    }


    return;
}

####
sub _process_cdhit_entries {
    my ($transcript_to_info_href, $transdecoder_acc_to_transcript_acc_href, $entries_aref, $cdhit_clusters_href, $cluster_counter) = @_;

    ## if ref entry is not PASA
    
    my @entries = @$entries_aref;
    my @transdecoder_structs;

    my $cluster_name = "cdhit-$cluster_counter";
    foreach my $entry (@entries) {
        my $transdecoder_acc = $entry->{acc};
        my $alignment_acc = $transdecoder_acc_to_transcript_acc_href->{$transdecoder_acc} or die "Error, no alignment obj found for transdecoder acc: $transdecoder_acc";
        my $alignment_struct = $transcript_to_info_href->{$alignment_acc} or die "Error, no struct for alignment $alignment_acc";

        my $transdecoder_orf_struct = $alignment_struct->{transdecoder_orfs}->{$transdecoder_acc} or die "Error, no transdecoder orf struct obtained for $transdecoder_acc";

        $transdecoder_orf_struct->{cdhit_clust_id} = $cluster_name;
        push (@transdecoder_structs, $transdecoder_orf_struct);
    }

    $cdhit_clusters_href->{$cluster_name} = [@transdecoder_structs];
        
}


sub populate_transdecoder_orfs {
    my ($transcript_to_info_href, $transdecoder_acc_to_transcript_acc_href, $trans_TD_pep_file) = @_;

    
    my $fasta_reader = new Fasta_reader($trans_TD_pep_file);
    while (my $seq_obj = $fasta_reader->next()) {

        my $acc = $seq_obj->get_accession();
        my $header = $seq_obj->get_header();
        my $pep_seq = $seq_obj->get_sequence();
        

        my $coding_length = length($pep_seq);

        my @header_components = split(/\s+/, $header);
        my $trans_info = pop @header_components;

        $trans_info =~ /^(\S+):(\d+)-(\d+)\(([\+\-])\)/ or confess "Error, cannot extract trans info from $trans_info";
        
        my $trans_acc = $1;
        my $lend = $2;
        my $rend = $3;
        my $orient = $4;
        
        if ($trans_acc =~ /Name=(S\d+)_(asmbl_\d+)/) {
            ## PASA transdecoder format
            ## definitely an unsightly hack
            $trans_acc = $1 . "-" . $2;
        }
        
        my $orf_struct = { lend => $lend,
                           rend => $rend,
                           orient => $orient,
                           aa_len => $coding_length,
                           acc => $acc,
                           pep => $pep_seq,
                       };

        my $transcript_struct = $transcript_to_info_href->{$trans_acc} or confess "Error, cannot get struct for $trans_acc";
        
        $transcript_struct->{transdecoder_orfs}->{$acc} = $orf_struct;
    
        $transdecoder_acc_to_transcript_acc_href->{$acc} = $trans_acc;
        
    }

    return;
}
    
sub populate_pasa_structure_info {
    my ($transcript_to_info_href, $pasa_gff3_file) = @_;
    
    my $alignment_indexer_href = {};
    
    my %scaffold_to_alignment_accs = &GFF3_alignment_utils::index_GFF3_alignment_objs($pasa_gff3_file, $alignment_indexer_href);

    foreach my $scaffold (keys %scaffold_to_alignment_accs) {
        
        my @align_accs = @{$scaffold_to_alignment_accs{$scaffold}};
        
        foreach my $align_acc (@align_accs) {
            my $alignment_obj = $alignment_indexer_href->{$align_acc};
            #print $alignment_obj->toString();
            

            my $acc = $alignment_obj->get_acc();
            my $orient = $alignment_obj->get_aligned_orientation();
            my $scaffold = $alignment_obj->{genome_acc};
            my ($lend, $rend) = sort {$a<=>$b} $alignment_obj->get_coords();
            my ($mlend, $mrend) = sort {$a<=>$b} $alignment_obj->get_mcoords();
            
            my $struct = { 
                
                is_transcript => 0,

                acc => $align_acc,
                length => $alignment_obj->{cdna_length},
                mapped => 1,
                valid => 1,
                alignment_structs => [  {  alignment => $alignment_obj->toToken(),
                                           valid => 1,
                                           scaffold => $scaffold,
                                           lend => $lend,
                                           rend => $rend,
                                           mlend => $mlend,
                                           mrend => $mrend,
                                       } ],
                
                transdecoder_orfs => {},
            };

            $transcript_to_info_href->{$align_acc} = $struct;
        }
    }

    return;
}
    
####
sub chimeric_transcript_breaks_orf {
    my ($struct, $transcript_to_info_href, $transdecoder_acc_to_transcript_acc) = @_;

    my @orfs = values %{$struct->{transdecoder_orfs}};
    
    my @alignments = sort {$a->{mlend}<=>$b->{mlend}} (@{$struct->{alignment_structs}});

    ## only two alignments in this mode.
    
    my $break_pos_left = $alignments[0]->{mrend};
    my $break_pos_right = $alignments[1]->{mlend};

    foreach my $orf (@orfs) {
        my ($orf_lend, $orf_rend) = ($orf->{lend}, $orf->{rend});
        
        if ($break_pos_left > $orf_lend && $break_pos_left < $orf_rend
            &&
            $break_pos_right > $orf_lend && $break_pos_right < $orf_rend) {

            
            return("chimeric break detected: trans_break($break_pos_left-$break_pos_right) orf($orf_lend-$orf_rend)");
            

        }
        
    }
    
    return("");
    
}

####
sub ORF_partially_aligned {
    my ($struct, $transcript_to_info_href, $transdecoder_acc_to_transcript_acc) = @_;

    my @orfs = values %{$struct->{transdecoder_orfs}};
    
    unless (@orfs) { 
        return("");
    }

    my $alignment = $struct->{alignment_structs}->[0];
    
    ## only two alignments in this mode.
    
    my $align_left_end = $alignment->{mlend};
    my $align_right_end = $alignment->{mrend};;
    
    ## just examining the longest ORF for now.
    @orfs = reverse sort {$a->{aa_len}<=>$b->{aa_len}} @orfs;
    
    my $longest_orf = shift @orfs;
    

    foreach my $orf ($longest_orf) { #@orfs) {
        my ($orf_lend, $orf_rend) = ($orf->{lend}, $orf->{rend});
        
        if (

            $orf_lend < $align_right_end && $orf_rend > $align_left_end
            &&

            ($orf_lend < $align_left_end-3 || $orf_rend > $align_right_end+3) # adjust for codon space

            ) {
            return("ORF only partially aligned (orf:$orf_lend-$orf_rend) vs. align:$align_left_end-$align_right_end");
            
        }
        
    }
    
    return("");
    
}


####
sub parse_cdna_span_from_alignment_text {
    my ($alignment) = @_;

    my @coords;
    while ($alignment =~ /\((\d+)\)/g) {
        push (@coords, $1);
    }

    @coords = sort {$a<=>$b} @coords;

    my $left_coord = shift @coords;
    my $right_coord = pop @coords;

    return($left_coord, $right_coord);
}





####
sub compare_to_PASA_cdhit_orfs {
    my ($struct, $transcript_to_info_href, $transdecoder_acc_to_transcript_acc_href, $cdhit_clusters_href) = @_;
    

    ## Does the transcript encode an ORF that is part of a CD-HIT cluster where it is the longer than any included PASA
    
    my %results = ( longer_than_pasa => 0,
                    shorter_than_pasa => 0,
                    same_as_pasa => 0 );
    
    foreach my $td_orf_acc (keys %{$struct->{transdecoder_orfs}}) {
        
        my $td_orf_struct = $struct->{transdecoder_orfs}->{$td_orf_acc};
        my $aa_len = $td_orf_struct->{aa_len};


        if (my $cdhit_cluster_id = $td_orf_struct->{cdhit_clust_id}) {
            
            my @clustered_tds = @{$cdhit_clusters_href->{$cdhit_cluster_id}};
            
            @clustered_tds = reverse sort {$a->{aa_len}<=>$b->{aa_len}} @clustered_tds;
            
            foreach my $clustered_td (@clustered_tds) {

                my $clustered_td_acc = $clustered_td->{acc};
                if ($clustered_td_acc eq $td_orf_acc) { next; }
        

        
                my $clustered_td_aa_len = $clustered_td->{aa_len};

                if ($clustered_td_acc =~ /PASA/) {
                    
                    my $transdecoder_trans_acc = $transdecoder_acc_to_transcript_acc_href->{$clustered_td_acc} or die "Error, no trans acc based on $clustered_td_acc";
                    my $transdecoder_trans_struct = $transcript_to_info_href->{$transdecoder_trans_acc} or die "Error, no trans struct based on $transdecoder_trans_acc";
                    
                    unless (&do_transcripts_have_overlapping_alignments($struct, $transdecoder_trans_struct)) { next; }
                                        
                    if ($clustered_td_aa_len == $aa_len) {
                        $results{same_as_pasa} = 1;
                    }
                    elsif ($clustered_td_aa_len > $aa_len) {
                        $results{shorter_than_pasa} = $aa_len/$clustered_td_aa_len * 100 unless ($results{shorter_than_pasa});
                        $results{shorter_than_pasa_note} = "aa_len: $aa_len vs. pasa_aa_len: $clustered_td_aa_len";
                        
                    }
                    elsif ($clustered_td_aa_len < $aa_len) {
                        $results{longer_than_pasa} = $aa_len/$clustered_td_aa_len * 100 unless ($results{longer_than_pasa});
                        $results{longer_than_pasa_note} = "aa_len: $aa_len vs. pasa_aa_len: $clustered_td_aa_len";
                    }
                }
            }
        }

                
    }
    
       
    return(%results);
    
    
}


####
sub PASA_SAME {
    my ($struct, $transcript_to_info_href, $transdecoder_acc_to_transcript_acc_href) = @_;
        

    my @pasa_accs_encapsulated = grep { $_->{acc} =~ /PASA|ANNOT/ } @{$struct->{cdhit_same_length}};

    my $alignment_struct = $struct->{alignment_structs}->[0];

    foreach my $pasa (@pasa_accs_encapsulated) {

        my $pasa_acc = $pasa->{acc};
        my $pasa_coding_len = $pasa->{aa_len};
        
        my $pasa_align_acc = $transdecoder_acc_to_transcript_acc_href->{$pasa_acc} or confess "Error, no transdecoder acc conversion for $pasa_acc";
        
        my $pasa_struct = $transcript_to_info_href->{$pasa_align_acc};

        my $pasa_alignment_struct = $pasa_struct->{alignment_structs}->[0];

        if (&do_alignments_overlap($alignment_struct, $pasa_alignment_struct)) {
            
            return("sufficiently represented by $pasa_acc w/ coding len: $pasa_coding_len");
        }
                
    }
    
       
    return("");
    
    
}

####
sub inferior_to_PASA {
    my ($struct, $transcript_to_info_href, $transdecoder_acc_to_transcript_acc_href) = @_;
    
    
    my @pasa_accs_encapsulating = grep { $_->{acc} =~ /PASA|ANNOT/ } @{$struct->{cdhit_encapsulated_by_longer}};

    my $alignment_struct = $struct->{alignment_structs}->[0];
    
    my $cdhit_max_aa_len = $alignment_struct->{cdhit_max_aa_len};

    foreach my $pasa (@pasa_accs_encapsulating) {

        my $pasa_acc = $pasa->{acc};
        my $pasa_coding_len = $pasa->{aa_len};
        
        my $pasa_align_acc = $transdecoder_acc_to_transcript_acc_href->{$pasa_acc} or confess "Error, no transdecoder acc conversion for $pasa_acc";
        
        my $pasa_struct = $transcript_to_info_href->{$pasa_align_acc};

        my $pasa_alignment_struct = $pasa_struct->{alignment_structs}->[0];

        if (&do_alignments_overlap($alignment_struct, $pasa_alignment_struct)) {
            
            return("already encapsulated by $pasa_acc w/ coding len: $pasa_coding_len");
        }
                
    }
    
       
    return("");
    
    
}


####
sub do_transcripts_have_overlapping_alignments {
    my ($structA, $structB) = @_;

    my @align_structs_A = @{$structA->{alignment_structs}};
    my @align_structs_B = @{$structB->{alignment_structs}};

    foreach my $alignA (@align_structs_A) {
        foreach my $alignB (@align_structs_B) {

            if (&do_alignments_overlap($alignA, $alignB)) {
                return(1);
            }
        }
    }

    return(0); # no overlaps

}



####
sub do_alignments_overlap {
    my ($align_struct_A, $align_struct_B) = @_;

    if ($align_struct_A->{scaffold} eq $align_struct_B->{scaffold}
        &&
        $align_struct_A->{lend} < $align_struct_B->{rend}
        &&
        $align_struct_A->{rend} > $align_struct_B->{lend} ) {

        return(1); # yes
    }
    else {
        return(0); # no
    }
}


####
sub summarize_orfs_found {
    my @orf_structs = @_;

    my $text = "";
    
    foreach my $orf_struct (@orf_structs) {
        if ($text) {
            $text .= "\t";
        }
        $text .= "ORF: " . $orf_struct->{lend} . "-" . $orf_struct->{rend} 
        . " len: " . $orf_struct->{aa_len} . "aa" 
            . ", Pep: " . $orf_struct->{pep};
    }
    
    return($text);
}

####
sub populate_chimera_jaccard_info {
    my ($transcript_to_info_href, $trans_chim_jaccard_summary_file) = @_;

    open (my $fh, $trans_chim_jaccard_summary_file) or die $!;

    while (<$fh>) {
        chomp;
        my ($acc, $trans_clip_pt, $jaccard, $single_count, $both_count, $map_left, $map_right) = split(/\t/);
        
        my $struct = $transcript_to_info_href->{$acc} or die "Error, no struct for $acc";
        
        $struct->{chimera_info} = { trans_clip_pt => $trans_clip_pt,
                                    jaccard => $jaccard,
                                    single_count => $single_count,
                                    both_count => $both_count,
                                    map_left => $map_left,
                                    map_right => $map_right,
                                };
        
    }

    close $fh;

    return;
}


####
sub populate_genome_seq_gaps {
    my ($gaps_file) = @_;

    my %scaff_to_gaps;

    open (my $fh, $gaps_file) or die $!;
    while (<$fh>) {
        chomp;
        my ($scaffold, $gap_lend, $gap_rend) = split(/\t/);

        my $gap_len = $gap_rend - $gap_lend + 1;
        if ($gap_len > $MIN_GAP_LEN) {
            push (@{$scaff_to_gaps{$scaffold}}, [$gap_lend, $gap_rend]);
        }
    }
    close $fh;

    return(%scaff_to_gaps);
}

####
sub near_sequencing_gap {
    my ($struct, $scaff_to_gap_positions_href) = @_;

    my @alignments = @{$struct->{alignment_structs}};
 

    my @gap_notes;
       
    foreach my $alignment (@alignments) {
        
        my $scaff = $alignment->{scaffold};
        my $lend = $alignment->{lend};
        my $rend = $alignment->{rend};
        
         if (my $gaps_list_aref = $scaff_to_gap_positions_href->{$scaff}) {
            
            foreach my $gap (@$gaps_list_aref) {

                my ($gap_lend, $gap_rend) = @$gap;
                
                if ($gap_lend < $rend && $gap_rend > $lend) {
                    push (@gap_notes, "seq gap:$scaff:$gap_lend-$gap_rend overlaps alignment:$scaff:$lend-$rend");
                }
                if ($lend > $gap_rend && $lend - $gap_rend <= $DIST_TO_GAP_REPORT) {
                    push (@gap_notes, "seq gap:$scaff:$gap_lend-$gap_rend within left range of alignment:$scaff:$lend-$rend");
                }
                if ($gap_lend > $rend && $gap_lend - $rend <= $DIST_TO_GAP_REPORT) {
                    push (@gap_notes, "seq gap:$scaff:$gap_lend-$gap_rend within right range of alignment:$scaff:$lend-$rend");
                }
            }
        }

    }

    if (@gap_notes) {
        return(join("\t", @gap_notes));
    }
    else {
        return("");
    }

}

####
sub populate_scaffold_lengths {
    my ($scaffold_lengths_file) = @_;

    my %scaff_lens;

    open (my $fh, $scaffold_lengths_file) or die $!;
    while (<$fh>) {
        chomp;
        my ($scaff, $len) = split(/\t/);
        $scaff_lens{$scaff} = $len;
    }
    close $fh;

    return(%scaff_lens);
}


####
sub near_scaffold_edge {
    my ($struct, $scaffold_lengths_href) = @_;

    my @alignments = @{$struct->{alignment_structs}};
    
    my @near_scaff_edge_notes;

    foreach my $alignment (@alignments) {
        my ($align_lend, $align_rend) = ($alignment->{lend}, $alignment->{rend});
        my ($align_mlend, $align_mrend) = ($alignment->{mlend}, $alignment->{mrend});
        my $scaff = $alignment->{scaffold};

        my $scaff_length = $scaffold_lengths_href->{$scaff} or die "Error, no length for scaff: $scaff";
        
        if ($align_lend <= $DIST_TO_SCAFF_EDGE_REPORT) {
            
            push (@near_scaff_edge_notes, "alignment:$scaff:$align_lend-$align_rend($align_mlend-$align_mrend) is adjacent to left edge of scaffold");
        }
        if ($scaff_length - $align_rend <= $DIST_TO_SCAFF_EDGE_REPORT) {
            
            push (@near_scaff_edge_notes, "alignment:$scaff:$align_lend-$align_rend($align_mlend-$align_mrend) is adjacent to right edge of scaffold: $scaff_length");
        }
    }
    
    return(@near_scaff_edge_notes);
    
}

####
sub has_overlapping_alignments {
    my ($struct) = @_;

    my @alignments = @{$struct->{alignment_structs}};

    unless (scalar @alignments > 1) { return (0); }

    for (my $i = 0; $i < $#alignments; $i++) {

        my $align_i = $alignments[$i];
        
        my $scaff_i = $align_i->{scaffold};
        my ($i_lend, $i_rend) = ($align_i->{lend}, $align_i->{rend});
        
        for (my $j = $i + 1; $j <= $#alignments; $j++) {

            my $align_j = $alignments[$j];
            
            my $scaff_j = $align_j->{scaffold};
            my ($j_lend, $j_rend) = ($align_j->{lend}, $align_j->{rend});

            if ($scaff_i eq $scaff_j
                &&
                $i_lend < $j_rend  && $i_rend > $j_lend) {

                return(1);
            }
        }
        

        
    }
    
    return(0); # no overlapping alignments.

}

####
sub populate_proteomics_info {
    my ($transcript_to_info_href, $proteomics_files) = @_;

    my @files = split(/,/, $proteomics_files);
    
    foreach my $file (@files) {

        open (my $fh, $file) or die "Error, cannot open file \"$file\"";
        while (<$fh>) {
            if (/^\#/) { next; } # header lines
            chomp;
            my ($acc, $plus_hits, $minus_hits) = split(/\t/);
            
            my $transcript_struct = $transcript_to_info_href->{$acc} or die "Error, no cdna struct for $acc";
            
            $transcript_struct->{proteomics_info}->{$file} = "+[$plus_hits];-[$minus_hits]";
            
        }
        close $fh;

    }

    return;
}
