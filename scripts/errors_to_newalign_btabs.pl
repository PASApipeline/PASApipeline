#!/usr/bin/env perl

use FindBin;
use lib ($FindBin::Bin);
use Pasa_init;
use strict;
use DBI;
use Getopt::Std;
use Mysql_connect;
use Ath1_cdnas;
use CdbTools;
use vars qw ($opt_r $opt_N $opt_M $opt_p $opt_s $opt_e $opt_n $opt_g $opt_d $opt_v $opt_h $opt_G $opt_b);
open (STDERR, "&>STDOUT");
&getopts ('M:p:s:e:n:g:dvhG:r:b:');
my $usage =  <<_EOH_;

Script locates the non-validating blat alignments and re-searches them using sim4.

############################# Options ###############################
# -M Mysql database/server ie. ("ath1_cdnas:haasbox")
# -p password file  (contains "username:password")
# -s cdna_seq_db (fasta file format)).
# -e prog containing error alignments [blat,sim4]
# -n next cdna alignment program to try. [sim4,geneseqer]
# -g geneseqer_org (only need if geneseqer is chosen as the -n option)
# -r geneseqer parameter file
# -d Debug
# -v verbose
# -h print this option menu and quit
# -G genomic_seq.db fasta file.
# -b bounds around failed alignment region to include in realignment (default 20kb)
#
###################### Process Args and Options #####################

_EOH_

    ;

if ($opt_h) {die $usage;}
my $MYSQLstring = $opt_M or die $usage;
my ($MYSQLdb, $MYSQLserver) = split (/:/, $MYSQLstring); 
my $genomic_seq_db = $opt_G or die $usage;
my $passwordinfo = $opt_p or die $usage;
my $DEBUG = $opt_d;
my $seqdb = $opt_s or die $usage;
my $delta = $opt_b || 20000;

my $error_prog = $opt_e or die $usage;
my $next_prog = $opt_n or die $usage;
my $geneseqer_org = $opt_g;
my $geneseqerprmfile = $opt_r;
if ($next_prog eq 'geneseqer' && (!$geneseqer_org || !$geneseqerprmfile))  { 
    die "Need to specify geneseqer_org and parameter file if using the geneseqer program";
}

our $SEE = $opt_v;
our $DB_SEE = $opt_v;

my %supported_error_progs = (sim4=>1, blat=>1, gmap=>1);
unless ($supported_error_progs{$error_prog}) { die "Sorry, $error_prog not currently supported under -e";}

my %supported_next_prog = (sim4=>1, geneseqer=>1);
unless ($supported_next_prog{$next_prog}) { die "Sorry, $next_prog not currently supported under -n";}


unless (-s $seqdb) {
    die "Sorry, can't locate $seqdb\n";
}

my ($user, $password) = split (/:/, $passwordinfo);

my ($dbproc) = &connect_to_db($MYSQLserver,$MYSQLdb,$user,$password);


my $completebtabfile = "$next_prog.${error_prog}_failed.btabs";
unlink ($completebtabfile) if -e $completebtabfile;
system "touch $completebtabfile"; #create empty file, in case there aren't any alignment so reprocess.  Next pipeline component expects this file.

## Query the failed searches:
my $query = "select c.annotdb_asmbl_id, cl.cdna_acc from cluster_link cl, cdna_link cdl, clusters c where c.cluster_id = cl.cluster_id and cl.cdna_acc = cdl.cdna_acc and cdl.prog = \"$error_prog\" and cdl.validate = 0 order by c.annotdb_asmbl_id";
my @results = &Mysql_connect::do_sql_2D ($dbproc, $query);
my $current_asmbl_id;
my $sequence = "";
foreach my $result (@results) {
    my ($asmbl_id, $acc) = @$result;
    print "Processing: $asmbl_id, $acc\n";
    
    my ($lend, $rend) = &get_failed_alignment_bounds($error_prog, $acc);
    

    if ($current_asmbl_id ne $asmbl_id) {
        $current_asmbl_id = $asmbl_id;
        $sequence = uc (&cdbyank_linear($asmbl_id, $genomic_seq_db));
    }
    
    my $tempseq = $acc;
    $tempseq =~ s/\W/_/g;
    $lend -= $delta;
    $lend = 0 if ($lend < 0);
    $rend += $delta; 
    ## pull out just the coorspanned region, minus a few kb on each side.
    my $subseq = substr ($sequence, $lend +1 -1, $rend - $lend + 1); #seq starts at nt downstream from lend; 
    
    $subseq =~ s/(\w{60})/$1\n/g; #convert to fasta format.
    
    open (GSEQ, ">gseq") or die;
    print GSEQ ">$asmbl_id\n$subseq\n";
    close GSEQ;
    
    my $cdna_seq = cdbyank_linear($acc, $seqdb);
    $cdna_seq =~ s/(\w{60})/$1\n/g; #convert to fasta format.
    open (CDNA, ">$tempseq") or die $!;
    print CDNA ">$acc\n$cdna_seq\n";
    close CDNA;
    
    
    my $btabfile = "$next_prog.btab";
    my $tmpoutput = "$next_prog.out";
    
    my $cmd = &prog_opts($tempseq, $next_prog) . " > $tmpoutput";
    print "$cmd\n" if $SEE;
    my $ret = system $cmd;
    die "$cmd\ndied" if $ret;
    
    
    my $cmd = "$ENV{PASAHOME}/scripts/${next_prog}_to_btab.pl <$tmpoutput >$btabfile";
    print "$cmd\n" if $SEE;
    my $ret = system $cmd;
    die "$cmd\ndied" if $ret;
    
    open (PROGFROM, $btabfile) or die;
    open (PROGTO, ">>$completebtabfile") or die;
    while (<PROGFROM>) {
        my @x = split (/\t/);
        $x[5] = $acc;
        $x[6] += $lend;
        $x[7] += $lend;
        my $out = join ("\t", @x);
        print PROGTO $out;
    }
    close PROGFROM;
    close PROGTO;
    exit if $DEBUG;
    unlink ($tempseq, $tmpoutput, $btabfile, "gseq");
}

sub prog_opts {
    my ($queryseq, $prog) = @_;
    if ($prog eq "sim4") {
        $prog = "sim4-mod"; #using the modified version of sim4 that's distributed with PASA, provided by Liliana Florea.  Only a slight output modification is all.
        
        
        return (" $prog $queryseq gseq A=4 ");
    }
    if ($prog =~ /geneseqer/i) {
        $prog = "GeneSeqer";
    }
    return (" $prog -s $geneseqer_org -E $queryseq -L gseq -m 1 -p $geneseqerprmfile "); #only reporting single best match.
}


####
sub get_failed_alignment_bounds {
    my ($prog, $cdna_acc) = @_;
    
    my $query = qq {  
        select min(a.lend), max(a.rend) from alignment a, cdna_link c 
            where c.align_id = a.align_id and c.cdna_acc = ? and c.prog = ? and c.validate = 0
        };
    
    my $result = &first_result_sql($dbproc, $query, $cdna_acc, $prog);
    return (@$result);
}



