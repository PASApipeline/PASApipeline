package main;
our $DEBUG;

package Eutilities_ncbi;

use strict;
use LWP::UserAgent;
use Data::Dumper;
use HTTP::Cookies;
use Carp qw (cluck longmess confess);

BEGIN {
   use Exporter;
   our (@ISA, @EXPORT);
   @ISA = qw (Exporter);
   @EXPORT = qw (egquery esearch esummary efetch efetch_batch elink
                 epost_uids epost_file link_post post_set print_summary 
                 print_links print_link_summaries );
   
   
   
}

my $delay = 0;
my $maxdelay = 3;

my $ua =  LWP::UserAgent->new;
$ua->cookie_jar(HTTP::Cookies->new);
$ua->timeout(1000);

#Contains the following subrotines:
#egquery
#esearch
#esummary
#efetch
#efetch_batch
#elink
#epost_uids
#epost_file
#link_post
#post_set
#print_summary
#print_links
#print_link_summaries

#*************************************************************

#************************************************************************

sub egquery {
# Performs EGQuery.
# Input: %params:
# $params{'term'} - Entrez query
# $params{'tool'} - tool name
# $params{'email'} - e-mail address
# Output = %results; keys are databases, values are UID counts

    my %params = @_;
    my $base = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/";
    my ($url, $raw);
    my @out;
    my $database;
    my %results;
    my ($begin, $end);
    
    sleep($delay);
    
    $url = $base . "egquery.fcgi?term=$params{'term'}";
    $url .= "&tool=$params{'tool'}&email=$params{'email'}";
    
    $begin = time;
    $raw = get($url);
    
    @out = split(/^/, $raw);
    
    foreach (@out) {
        
        if (/<DbName>(.*)<\/DbName>/) { $database = $1; }
        if (/<Count>(\d+)<\/Count>/) { $results{$database} = $1; }
        
    }
    
    $end = time;
    $delay = $maxdelay - ($end - $begin);
    if ($delay < 0) { $delay = 0; }
    
    return(%results);
    
}

#*********************************************************************

sub esearch {
# Performs ESearch.
# Input: %params
# $params{'db'} - database
# $params{'term'} - Entrez query
# $params{'usehistory'} (y/n) - flag for using the Entrez history server
# $params{'retstart'} - first item in results list to display (default = 0)
# $params{'retmax'} - number of items in results list to display (default = 20)
# $params{'WebEnv'} - Web Environment for accessing existing data sets
# $params{'tool'} - tool name
# $params{'email'} - e-mail address
#
# Output: %results: keys are 'count', 'query_key', 'WebEnv', 'uids'
# $results{'uids'} is an array
   
    my %params = @_;
    my $base = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/";
    my ($url, $raw);
    my @out;
    my %results;
    my ($begin, $end);
    
    if ($DEBUG) {
        print "ESEARCH: " . Dumper (\%params);
    }
    
    sleep($delay);
    
    $url = $base . "esearch.fcgi?db=$params{'db'}&term=$params{'term'}";
    $url .= "&usehistory=$params{'usehistory'}&WebEnv=$params{'WebEnv'}";
    $url .= "&retstart=$params{'retstart'}&retmax=$params{'retmax'}";
    $url .= "&tool=$params{'tool'}&email=$params{'email'}";
    
    $begin = time;
    $raw = get($url);
    
    print "RAW: $raw\n" if $DEBUG;
    
    $raw =~ /<Count>(\d+)<\/Count>/s;
    $results{'count'} = $1;
    $raw =~ /<QueryKey>(\d+)<\/QueryKey>.*<WebEnv>(\S+)<\/WebEnv>/s;
    $results{'query_key'} = $1 if ($params{'usehistory'} eq 'y');
    $results{'WebEnv'} = $2;
    @out = split(/^/, $raw);
    
    foreach (@out) {
        if (/<Id>(\d+)<\/Id>/) { push (@{$results{'uids'}}, $1); }
    }
    
    $end = time;
    $delay = $maxdelay - ($end - $begin);
    if ($delay < 0) { $delay = 0; }
    
    return(%results);
    
}

#****************************************************************

sub esummary {
    # Performs ESummary. 
    # Input: %params:
    # $params{'db'} - database
    # $params{'id'} - UID list 
    # $params{'query_key'} - query_key
    # $params{'WebEnv'} - web environment
    # $params{'retstart'} - first DocSum to retrieve
    # $params{'retmax'} - number of DocSums to retrieve
    # $params{'tool'} - tool name
    # $params{'email'} - e-mail address 
    # 
    # Output: %results: $results{id}{item} = value
    # where id = UID, item = Item Name
    
    my %params = @_;
    my $base = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/";
    my ($url, $raw);
    my @out;
    my $id;
    my %results;
    my ($begin, $end);
    
    sleep($delay);
    
    $url = $base . "esummary.fcgi?db=$params{'db'}&retstart=$params{'retstart'}";
    $url .= "&retmax=$params{'retmax'}";
    
    if ($params{'query_key'}) {
        $url .= "&query_key=$params{'query_key'}&WebEnv=$params{'WebEnv'}";   
    }
    else {
        $url .= "&id=$params{'id'}";
    }
    
    $url .= "&tool=$params{'tool'}&email=$params{'email'}";
    
    $begin = time;
    
    # print "URL: $url\n";
    
    $raw = get($url);
    
    # print "RAW: $raw\n";
    
    @out = split(/^/, $raw);
    
    foreach (@out) {
	
        $id = $1 if (/<Id>(\d+)<\/Id>/);
        if (/<Item Name="(.+)" Type=.*>(.+)<\/Item>/) {
            $results{$id}{$1} = $2;
        }
        
    }
    
    $end = time;
    $delay = $maxdelay - ($end - $begin);
    if ($delay < 0) { $delay = 0; }
    
    return(%results);
    
}

#****************************************************************

sub efetch {
    # Performs EFetch. 
    # Input: %params:
    # $params{'db'} - database
    # $params{'id'} - UID list 
    # $params{'query_key'} - query key
    # $params{'WebEnv'} - web environment
    # $params{'retmode'} - output data format
    # $params{'rettype'} - output data record type
    # $params{'retstart'} - first record in set to retrieve
    # $params{'retmax'} - number of records to retrieve
    # $params{'seq_start'} - retrieve sequence starting at this position
    # $params{'seq_stop'} - retrieve sequence until this position
    # $params{'strand'} - which DNA strand to retrieve (1=plus, 2=minus)
    # $params{'complexity'} - determines what data object to retrieve
    # $params{'tool'} - tool name
    # $params{'email'} - e-mail address
    #
    # Output: $raw; raw EFetch output

    my %params = @_;
    my $base = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/";
    my ($url, $raw);
    my ($begin, $end);
    
    sleep($delay);
    
    
    print STDERR Dumper (\%params) if $DEBUG;
    

    $url = $base . "efetch.fcgi?db=$params{'db'}";
    
    if ($params{'query_key'}) {
        $url .= "&query_key=$params{'query_key'}&WebEnv=$params{'WebEnv'}";
    }
    else {
        $url .= "&id=$params{'id'}";
    }
    
    $url .= "&retmode=$params{'retmode'}&rettype=$params{'rettype'}";
    $url .= "&retstart=$params{'retstart'}&retmax=$params{'retmax'}";
    $url .= "&seq_start=$params{'seq_start'}&seq_stop=$params{'seq_stop'}";
    $url .= "&strand=$params{'strand'}&complexity=$params{'complexity'}";
    $url .= "&tool=$params{'tool'}&email=$params{'email'}";
    
    $raw = "Error:";
    
    print STDERR "URL: $url\n" if $DEBUG;
    
    while ($raw =~ /Error\:/) {
        
        $begin = time;
        
        $raw = get($url);
        #print "Getting $url\n";
        
        $end = time;
        $delay = $maxdelay - ($end - $begin);
        if ($delay < 0) { $delay = 0; }
        
        last; # breaking the error checking loop for now.

        if ($raw !~ /Error\:/) {
            
            last; # got what we came for
        }
        else {
            print STDERR $raw;
            sleep ($delay); #wait a short time before trying again
        }
    }
    
    if ($DEBUG) {
        open (my $fh, ">retstart_" . $params{retmax} . "_" . $params{retstart}) or die $!;
        print $fh $raw;
        close $fh;
    }
    
    return($raw);
    
}

#****************************************************************

sub efetch_batch {
    # Uses efetch to download a large data set in 500 record batches
    # The data set must be stored on the History server
    # Input: %params:
    # $params{'db'} - link to database
    # $params{'query_key'} - query key
    # $params{'WebEnv'} - web environment
    # $params{'retmode'} - output data format
    # $params{'rettype'} - output data record type
    # $params{'seq_start'} - retrieve sequence starting at this position
    # $params{'seq_stop'} - retrieve sequence until this position
    # $params{'strand'} - which DNA strand to retrieve (1=plus, 2=minus)
    # $params{'complexity'} - determines what data object to retrieve
    # $params{'tool'} - tool name
    # $params{'email'} - e-mail address
    # $params{'retmax'} -number of entries fetched per ncbi eutilities call (default 500)
    # $params{'filehandle'} filehandle to write data to, or uses stdout
    # $params{'returnOutput'}  set to 1 if you want the results returned as a string
    #
    # Output: none (default) or all the output if 'returnOutput' flag is turned on.
    
    my %params = @_;
    my $base = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/";
    my ($url, $raw);
    my ($begin, $end);
    my %results;
    my ($count, $limit);
    my ($retstart, $retmax);
    
    my $filehandle = $params{filehandle};
    
    $retmax = $params{retmax} || 500;
    
    #first use ESearch to determine the size of the dataset
    
    $params{'term'} = "%23" . "$params{'query_key'}";
    $params{'usehistory'} = 'y';
    
    %results = esearch(%params);
    
    $count = $results{'count'};
    $params{'retmax'} = $retmax;
    
    print STDERR "Retrieving $count records...\n";
    
    for ($retstart = 0; $retstart < $count; $retstart += $retmax) {
        
        sleep($delay);
        $params{'retstart'} = $retstart;
        $begin = time;
        my $fetch_results = efetch(%params);
        print $filehandle $fetch_results if $filehandle;
        
        if ($params{returnOutput}) {
            $raw .= $fetch_results;
        }
        
        if ($retstart + $retmax > $count) { $limit = $count; }
        else { $limit = $retstart + $retmax; }
        
        print STDERR "Received records $retstart - $limit.\n";
        $end = time;
        $delay = $maxdelay - ($end - $begin);
        if ($delay < 0) { $delay = 0; }
    }
    
    return ($raw);
    
}

#****************************************************************

sub elink {
    # Performs ELink.
    # Input: %params:
    # $params{'dbfrom'} - link from database
    # $params{'db'} - link to database
    # $params{'id'} - array of UID lists
    # $params{'query_key'} - query key
    # $params{'WebEnv'} - web environment)
    # $params{'term'} - Entrez term used to limit link results
    # $params{'tool'} - tool name
    # $params{'email'} - e-mail address
    #
    # Output: %links:
    # @{$links{'from'}{$set}} = array of input UIDs in set $set
    # @{$links{'to'}{$db}{$set}} = array of linked UIDs in $db in set $set
    # where $set = integer corresponding to one &id parameter
    # value in the ELink URL

    my %params = @_;
    my $base = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/";
    my ($url, $raw);
    my ($line, $getdata, $getid, $link, $id, $set);
    my @out;
    my @link_ids;
    my $ids;
    my %results;
    my $db;
    my ($begin, $end);

    sleep($delay);

    $set = 0;
    
    $url = $base . "elink.fcgi?dbfrom=$params{'dbfrom'}&db=$params{'db'}";
    $url .= "&term=$params{'term'}";
    
    if ($params{'query_key'}) {
        
        $url .= "&query_key=$params{'query_key'}&WebEnv=$params{'WebEnv'}";
        
    }
    else {
        
        foreach $ids (@{$params{'id'}}) {
            $url .= "&id=$ids";
        }
    }
    $url .= "&tool=$params{'tool'}&email=$params{'email'}";
    
    $begin = time;
    $raw = get($url);
    
    print "link output raw:\n$raw\n" if $DEBUG;
    
    @out = split(/^/,$raw);
    
    $getdata = 0;
    
    foreach $line (@out) {
        
        #check for input UIDs
        $getid = 1 if ($line =~ /<IdList>/);
        if ($getid) {
            
            push (@{$results{'from'}{$set}}, $1) if ($line =~ /<Id>(\d+)<\/Id>/);
            
        }
        $getid = 0 if ($line =~ /<\/IdList>/);
        
        #check for linked UIDs   
        if ($line =~ /<DbTo>(\S+)<\/DbTo>/) {
            $db = $1;
            $getdata = 1;
        }
        $getdata = 0 if ($line =~ /<\/LinkSetDb>/);
        
        if ($getdata) {
            push (@{$results{'to'}{$db}{$set}}, $1) if ($line =~ /<Id>(\d+)<\/Id>/);
        }
        
        if ($line =~ /<\/LinkSet>/) {
            $getdata = 0;
            $set++;
        }
        
    }
    
    $end = time;
    $delay = $maxdelay - ($end - $begin);
    if ($delay < 0) { $delay = 0; }
    
    return(%results);
    
}

#*********************************************************************

sub epost_uids {
    # Performs EPost, placing UIDs in the URL.
    # Input: %params:
    # $params{'db'} - database
    # $params{'id'} - list of UIDs
    # $params{'tool'} - tool name
    # $params{'email'} - e-mail address
    #
    #Output: %results: keys are 'WebEnv' and 'query_key'

    my %params = @_;
    my $base = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/";
    my ($url, $raw);
    my ($begin, $end);
    my %results;
    sleep($delay);
    
    $url = $base . "epost.fcgi?db=$params{'db'}&id=$params{'id'}";
    $url .= "&tool=$params{'tool'}&email=$params{'email'}";
    
    $begin = time;
    $raw = get($url);
    
    $raw =~ /<QueryKey>(\d+)<\/QueryKey>.*<WebEnv>(\S+)<\/WebEnv>/s;
    $results{'query_key'} = $1;
    $results{'WebEnv'} = $2;
    
    $end = time;
    $delay = $maxdelay - ($end - $begin);
    if ($delay < 0) { $delay = 0; }
    
    return(%results);
    
}

#*********************************************************************

sub epost_file {
    # Performs EPost, accepts input from file. 
    # Input file must have one UID per line.
    # Input: %params:
    # $params{'db'} - database
    # $params{'id'} - filename containing a list of UIDs
    # $params{'tool'} - tool name
    # $params{'email'} - e-mail address
    #
    # Output: %results: keys are 'WebEnv' and 'query_key'

    my %params = @_;
    my $base = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/";
    my $uids;
    my @list;
    my ($begin, $end);
    my %results;
    
    sleep($delay);
    
    #read input file of UIDs, one per line
    my $file = $params{'id'};
        
    my @list = `cat $file`;
    chomp @list;
    unless (@list) {
        die "Errror, no uids in the list of file $file\n";
    }
    
    $params{'id'} = join (',', @list);
    
    unless ($params{id}) {
        die "Error, no IDs to process ";
    }
    
    %results = post_set(%params);
    
    $end = time;
    $delay = $maxdelay - ($end - $begin);
    if ($delay < 0) { $delay = 0; }
    
    return(%results);
    
}

#************************************************************

sub link_post {
    
    # Uses EPost to post ELink results on the History server
    # Input: %results output of sub elink
    # Output: %results: $results{$set}{$db}{'query_key'}
    #		    $results{$set}{$db}{'WebEnv'}
    
    my %linkout = @_;
    my $base = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/";
    my ($url, $raw, $uids, $key, $db);
    my (%params, %postout, %results);
    
    
    foreach $db (keys %{$linkout{'to'}}) {
        foreach $key (sort keys %{$linkout{'from'}}) {
            
            $params{'db'} = $db;
            $uids = join (',', @{$linkout{'to'}{$db}{$key}} );
            
            $params{'id'} = $uids;
            %postout = post_set(%params);
            $results{$key}{$db}{'query_key'} = $postout{'query_key'};
            $results{$key}{$db}{'WebEnv'} = $postout{'WebEnv'};
            
        }
        
    }
    
    return (%results);
    
}

#***********************************************************

sub post_set {
    
    my (%params) = @_;
    my $base = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/";
    my ($url_params, $raw, $url);
    my %results;
    
    $url_params = "db=$params{'db'}&id=$params{'id'}";
    $url_params .= "&tool=$params{'tool'}&email=$params{'email'}";
    
    $url = $base . "epost.fcgi";
    
    #create user agent
    my $ua = new LWP::UserAgent;
    $ua->agent("epost_file/1.0 " . $ua->agent);
    
    #create HTTP request object
    my $req = new HTTP::Request POST => "$url";
    $req->content_type('application/x-www-form-urlencoded');
    
    print "post_set() URL_params: $url_params\n" if $DEBUG;
    
    $req->content("$url_params");
    
    
    #post the HTTP request
    $raw = $ua->request($req);
    
    $raw->content =~ /<QueryKey>(\d+)<\/QueryKey>.*<WebEnv>(\S+)<\/WebEnv>/s;
    $results{'query_key'} = $1;
    $results{'WebEnv'} = $2;
    
    return (%results);
    
}

#***********************************************************

sub print_summary {
    
    # Input: %results output from sub esummary
    
    my %results = @_;
    my $id;
    
    foreach $id (sort keys %results) {
        
        print "\nID $id:\n" if $DEBUG;
        foreach (sort keys %{$results{$id}}) {
            print "$_: $results{$id}{$_}\n" if $DEBUG;
        }
    }
    
}

#***********************************************************

sub print_links {
    
    # Input: %results output from sub elink
    
    my %results = @_;
    my ($key, $db);
    
    foreach $key (sort keys %{$results{'from'}}) {
        print "Links from: " if $DEBUG;
        foreach (@{$results{'from'}{$key}}) {
            print "$_ " if $DEBUG;
        }
        foreach $db (keys %{$results{'to'}}) {
            print "\nto $db:" if $DEBUG;
            foreach (@{$results{'to'}{$db}{$key}}) {
                print "$_ " if $DEBUG;
            }
        }
        print "\n***\n" if $DEBUG;
    }
    
}

#**********************************************************

sub print_link_summaries {
    
    # Input: %results output from sub link_post
    
    my %results = @_;
    my (%params,%docsums);
    my ($db, $key);
    
    foreach $key ( sort keys %results ) {
        
        print "Links from set $key\n" if $DEBUG;
        foreach $db (keys %{$results{$key}} ) {
            
            $params{'db'} = $db;
            $params{'WebEnv'} = $results{$key}{$db}{'WebEnv'};
            $params{'query_key'} = $results{$key}{$db}{'query_key'};
            %docsums = esummary(%params);
            print "$db\n\n" if $DEBUG;
            print_summary(%docsums);
            print "\n" if $DEBUG;
        }
    }
    
}


sub get {
    my $url = shift;

    my $content = "";
    
    while (1) {
        my $response = $ua->get($url);
        if ($response->is_success) {
            $content = $response->content;
            last;
        }
        else {
            print STDERR $response->status_line;
            print STDERR "\t-trying again.\n";
        }
    }
    
    return ($content);
    
}

1; #EOM


