#!/usr/bin/env perl

use strict;
use CGI;

my $html_flag = (@ARGV) ? 1:0;

my $start_html = 0;

my $URL = "http://wild/tigr-scripts/euk_manatee/shared/ORF_infopage.cgi?db=ath1&orf";

my $cgi;
if ($html_flag) {
    $cgi = new CGI();
}

$/ = "\n>";
while (<STDIN>) {
    my $error;
    my $fastaSeq = $_;
    $fastaSeq =~ s/>//g;
    split (/\n/, $fastaSeq);
    #print $fastaSeq;
    my $header = shift @_ ;
    my $sequence = join ("", @_);
    #print "$sequence\n\n";
    my $numstops = 0;
    while ($sequence =~ /\*/g) {
	$numstops++;
    }
    if ($numstops >1 || !($sequence)) {
	$error .= "\t\*$numstops";
    } 
    unless ($sequence =~ /^m/i) {
	$error .= "\tDoesn't start with M.";
    }
    unless ($sequence =~ /\*$/) {
	$error .= "\tNo stop codon.";
    }
    if ($error) {
	if ($html_flag) {
	    if (!$start_html) {
		$start_html = 1;
		print $cgi->start_html();
		print "<table border=1>\n";
	    }
	    my @stuff = split (/\s+/, $header);
	    my $acc = shift @stuff;
	    my $url = "<a href=\"$URL=$acc\" target=_blank\">$acc</a>";
	    $header =~ s/$acc/$url/;
	    print "<tr><td>$header<br><pre>ERRORS:\n$error</pre></td></tr>\n";
	} else {
	    print "$header\tERRORS: $error\n\n";
	}
    }
}


if ($html_flag) {
    print "</table>" . $cgi->end_html();
}

