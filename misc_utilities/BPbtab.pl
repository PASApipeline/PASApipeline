#!/usr/bin/env perl


=head1 NAME
    
    BPbtab

    script parses WU-BLAST or NCBI-BLAST output files into BTAB format where each HSP is reported as a single line with tab-delimited fields.

    The original BTAB program is described here:

    Dubnick M. (1992) Btab--a Blast output parser. Comput Appl Biosci 8(6):601-2

    
    A Perl version which emulates the functionality of btab is provided here by BPbtab.  BPbtab relies exclusively on the BioPerl 1.4 Bio::SearchIO functionality.


    
=cut


=head1 USAGE

    Standard input is parsed and written to standard output.

    BPbtab < blast.output > blast.output.btab

=cut


use strict;
use Bio::SearchIO;

my $in = new Bio::SearchIO(-format => 'blast', 
			   -fh   => \*STDIN);


# parse each blast record:
while( my $result = $in->next_result ) {
    
    # parse each hit per record.
    while( my $hit = $result->next_hit ) {
        
        # a hit consists of one or more HSPs
        while( my $hsp = $hit->next_hsp ) {
     
            eval {
                &process_blast_match ($result, $hit, $hsp);
            };
            if ($@) {
                print STDERR "Error, couldn't process blast entry.\n";
            }
        }
    }
}


exit(0);


####
sub process_blast_match {
    my ($result, $hit, $hsp) = @_;
    
    my @x;
    $x[0] = $result->query_name();
    # date
    $x[2] = $result->query_length();
    $x[3] = $hsp->algorithm();
    $x[4] = $result->database_name();
    $x[5] = $hit->name();
    $x[6] = $hsp->start('query');
    $x[7] = $hsp->end('query');
    my $queryStrand = $hsp->strand('query');
    if ($queryStrand == -1) {
        ($x[6], $x[7]) = ($x[7], $x[6]);
    }
    
    $x[8] = $hsp->start('hit');
    $x[9] = $hsp->end('hit');
    my $hitStrand = $hsp->strand('hit');
    if ($hitStrand == -1) {
        ($x[8], $x[9]) = ($x[9], $x[8]);
    }
    
    $x[10] = sprintf ("%.1f", $hsp->percent_identity());   
    
    my $similarity = $hsp->frac_conserved('total') * 100; 
    $x[11] = sprintf("%.1f", $similarity);
    $x[12] = $hsp->score();
    $x[13] = $hsp->bits();
    
    $x[15] = $hit->description();
    
    $x[16] = ( ($hsp->query->frame + 1) * $hsp->query->strand); #blast frame (1, 2, 3, -1, -2, -3).
    
    
    my $strandDescript = "null";
    if ($queryStrand == 1) {
        $strandDescript = "Plus";
    } elsif ($queryStrand == -1) {
        $strandDescript = "Minus";
    }
    $x[17] = $strandDescript;
    
    $x[18] = $hit->length();
    $x[19] = $hsp->evalue();
    $x[20] = $hsp->pvalue();
    
    my $outline = join ("\t", @x);
    print "$outline\n";
    
}


=head1 DESCRIPTION

                                   Fields


The parsed BLAST output is presented as single line HSP descriptions, with tab-delimited fields in the following order:


[0]      Query Sequence Name

[2]      Query Sequence Length

[3]      Search Method  --  Blast family application name

[4]      Database Name

[5]      Subject Sequence Name  --  Database entry name

[6],[7]   Query Left End, Query Right End  --  The endpoints of the part
               of the query sequence which Blast aligns with the subject sequence.

[8],[9]  Subject Left End, Subject Right End  --  The endpoints of the part 
               of the subject sequence which Blast aligns with the query sequence.

[10]     Percent Identity  --  The fraction of residues which are absolute
            matches between the query and subject sequence, expressed in
            percent.  

[11]     Percent Similarity  --  The fraction of residues which are exact or 
            similar matches between the query and subject sequence, expressed in
            percent. 

[12]    HSP score    

[13]    Bits score


[15]    Description  --  A freeform text field which contains the biological
        description field from the database for the subject sequence.  If
        this text occupies more than one line in the Blast output file, the
        NewLines are replaced by spaces.  Commas may occur in this field even
        if they are the field separator character, because this is the last
        field in the record.

[16]    Query Frame (1, 2, 3, -1, -2, -3)

[17]    Query Strand  --  Plus, Minus or null

[18]    DB sequence length

[19]    Expect -- expected value

[20]    P-Value  --  Poisson ratio


** Note ** Intervening field positions which are not described are not currently supported.  These remain to support compatibility with other existing btab implementations.

=cut

