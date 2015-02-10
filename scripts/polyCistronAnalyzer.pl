#!/usr/bin/env perl

use Mysql_connect;
use DBI;
use Getopt::Long qw(:config no_ignore_case bundling);
use Longest_orf;
use Data::Dumper;
use Fasta_reader;
use strict;
require "overlapping_nucs.ph";


use vars qw ($MIN_PROTEIN_LENGTH $MAX_OVERLAP_PERCENT $opt_h $opt_P $opt_T $opt_M $opt_p);

&GetOptions ( "MIN_CDS_LENGTH=s" => \$MIN_PROTEIN_LENGTH,
	      "MAX_OVERLAP_PERCENT=s" => \$MAX_OVERLAP_PERCENT,
	      "h" => \$opt_h,
	      "P=s" => \$opt_P,
	      "T=s" => \$opt_T,
	      "M=s" => \$opt_M,
	      "p=s", \$opt_p);


my $usage = <<_EOH;


Finds examples of likely polyCistronic mRNAs by looking for those cDNAs that have more than one 
long ORF with minimal overlap

Also, identifies potential PASA chimeric assemblies defined as those which appear polycistronic but lack any 
composite potentially polycistronic mRNA.


###########################################################################################
## --MIN_PROTEIN_LENGTH   :minimum size of a protein encoded within the cDNA (default 300)
## --MAX_OVERLAP_PERCENT :two overlapping CDSs are allowed to overlap by this max percent (default 10)
##
## -P PASA assemblies (fasta file)
## -T transcripts (fasta file)
## -M Mysql database/server ie. ("ath1_cdnas:haasbox")
## -p passwordinfo  (contains "username:password")
##
## -h help menu
###########################################################################################


_EOH

    ;


if ($opt_h) {
    die $usage;
}
my $MYSQLstring = $opt_M or die $usage;
my ($MYSQLdb, $MYSQLserver) = split (/:/, $MYSQLstring); 
my $passwordinfo = $opt_p or die $usage;
my ($user, $password) = split (/:/, $passwordinfo);
my ($dbproc) = &Mysql_connect::connect_to_db($MYSQLserver,$MYSQLdb,$user,$password);



my $MIN_CDS_LENGTH = (3*$MIN_PROTEIN_LENGTH) ||  900;

unless (defined $MAX_OVERLAP_PERCENT) {
    $MAX_OVERLAP_PERCENT = 10;
}


main: {

    my %PASA_polycistrons = &find_polyCistrons($opt_P);
        
    my %Transcript_polycistrons = &find_polyCistrons($opt_T);
    
    my %asmbl_link = &get_asmbl_link();
    
    foreach my $pasa_acc (keys %PASA_polycistrons) {
	
	my $transcript_supported = 0;

	my $cdna_list_aref = $asmbl_link{$pasa_acc};
	foreach my $cdna_acc (@$cdna_list_aref) {
	    if ($Transcript_polycistrons{$cdna_acc}) {
		$transcript_supported = 1;
		last;
	    }
	}

	if ($transcript_supported) {
	    print "$pasa_acc\tPolyCistron\n";
	} else {
	    print "$pasa_acc\tPasaChimera\n";
	}
    }
        
    exit(0);
}



####
sub find_polyCistrons {
    my $fastaFile = shift;
    my $reader = new Fasta_reader($fastaFile);
 
    my %polyCistrons;

    my $polyCistron_outfile = "$fastaFile.polyCistrons";
    if (-s $polyCistron_outfile) {
	## already ran, just read the results
	print STDERR "Already have file $polyCistron_outfile, reading it now.\n";
	open (POLY, $polyCistron_outfile) or die "Cannot open $polyCistron_outfile";
	while (<POLY>) {
	    if (/^\/\//) {
		my @x = split (/\s+/);
		my $acc = $x[1];
		$polyCistrons{$acc} = 1;
	    }
	}
	close POLY;
    } 

    else {

	open (POLY, ">$polyCistron_outfile") or die "Cannot open file $polyCistron_outfile\n";
	
	while (my $seqEntry = $reader->next()) {
	    my $sequence = $seqEntry->get_sequence();
	    my $accession = $seqEntry->get_accession();
	    my $header = $seqEntry->get_header();
	    
	    my $longest_Orf = Longest_orf->get_longest_orf($sequence);
	    
	    my @orfs = $longest_Orf->orfs();
	    if ($#orfs > 0) { #more than one ORF found:
		my $outtext = "";
		my $first_orf = shift @orfs;
		my $second_orf = shift @orfs;
		my $f_length = $first_orf->{length};
		my $s_length = $second_orf->{length};
		if ($f_length < $MIN_CDS_LENGTH || $s_length < $MIN_CDS_LENGTH) {
		    next;
		}
		
		## See if they overlap:
		my ($f_lend, $f_rend) = sort {$a<=>$b} ($first_orf->{start}, $first_orf->{stop});
		my ($s_lend, $s_rend) = sort {$a<=>$b} ($second_orf->{start}, $second_orf->{stop});
		$outtext .= "// $header\n";
		if ($f_lend < $s_rend && $f_rend > $s_lend) { #overlap
		    $outtext .=  "Overlapping ORFs:\t";
		    my $overlap_length = nucs_in_common($f_lend, $f_rend, $s_lend, $s_rend);
		    my $f_percent_overlap = $overlap_length/$f_length*100;
		    my $s_percent_overlap = $overlap_length/$s_length*100;
		    if ($f_percent_overlap > $MAX_OVERLAP_PERCENT || $s_percent_overlap > $MAX_OVERLAP_PERCENT) {
			next;
		    }
		    
		} else {
		    $outtext .=  "Non-overlapping ORFs:\t";
		}
		$outtext .= " ($f_lend, $f_rend) vs. ($s_lend, $s_rend)\n";
		
		print POLY "$outtext";
		$polyCistrons{$accession} = 1;
	    }
	    
	}
	close POLY;
    }
    return (%polyCistrons);
}


####
sub get_asmbl_link {
    my $query = "select asmbl_acc, cdna_acc from asmbl_link";
    my @results = &do_sql_2D($dbproc, $query);
    
    my %asmbl_link;
    foreach my $result (@results) {
	my ($asmbl_acc, $cdna_acc) = @$result;
	my $list_ref = $asmbl_link{$asmbl_acc};
	unless (ref $list_ref) {
	    $list_ref = $asmbl_link{$asmbl_acc} = [];
	}
	push (@$list_ref, $cdna_acc);
    }

    return (%asmbl_link);
}
