#!/usr/local/bin/perl
#
# Mailer.pm
#
# Mailer - basic module for sending a simple (non-encoded) e-mail through a SMTP server
# for internal usage- so the default parameters (e.g. SMTP server) are set 
# to TIGR defaults
#
# It doesn't need any external shell program (like mail) - instead it connects to a mail server
# directly using sockets from Perl.
#
# Exports only one function: send_mail. Usage:
#
# send_mail {[to=>'<dest>', subj=>'<subj>', body=>'<body>']};
# optional hash elements: 
#      from=>'<addr>'  - default is <current_user>@<current_host>
#      smtp=>'<smtp>'  - default is: mail.<domain>
#      file=>'<file>'  - if this is given, body entry may be missing, but if body 
#              is given too, the content of <file> is just sent after the 'body'

use FileHandle;
use Socket;
use strict 'vars';
#default values (change them for your environment)
my $def_domain = `grep "^domain" /etc/resolv.conf`;
if ($def_domain) {
  chomp($def_domain);
  ($def_domain)=($def_domain =~ m/^\S+\s+(\S+)/);
 }
else { #fallback on some hardcoded value
 $def_domain = 'tigr.org';
 } 
my $def_host=`hostname`;
chomp($def_host);
$def_host.='.'.$def_domain unless $def_host =~ m/\.$def_domain/;
my $def_smtp='mail.'.$def_domain; #change this to your SMTP server
return 1;

#------------------------------------------------------------------------------------------

sub send_mail {
 return "Incorrect call to send_mail: a hash is expected." if (ref $_[0] ne 'HASH');
 my $hash=shift;
 $hash->{'smtp'}=$def_smtp unless defined $hash->{'smtp'};
 my $file;
 if (defined($hash->{file})) {
   #warning: assumes it's a decent, short text file!
   local $/=undef; #read whole file
   open(ADDFILE, '<'.$hash->{file}) || return "Cannot open file ".$hash->{file}."\n";
   $file=<ADDFILE>;
   close ADDFILE;
   }
 my $proto=(getprotobyname('tcp'))[2];
 my $port = getservbyname('smtp', 'tcp');
 #my $boundary = 'Message-Boundary-19990614';
 $hash->{'from'}=$ENV{'USER'}.'@'.$def_host unless defined($hash->{'from'});
 my $to=$hash->{'to'};
 unless ($to) {
    $hash->{'to'}=$ENV{'USER'}.'@'.$def_domain;
    }
   else {
    $hash->{'to'}=$to.'@'.$def_domain unless $to =~ m/@/;
    }
 $hash->{'cc'} =~ s/\s+/ /g if (defined $hash->{'cc'}); # pack spaces and add comma
 $hash->{'cc'} =~ s/,,/,/g if (defined $hash->{'cc'});
 $hash->{'smtp'} =~ s/^\s+//g; # remove spaces around $smtp
 $hash->{'smtp'} =~ s/\s+$//g;
 my $client=defined($hash->{'client'})?$hash->{'client'}:$def_host;
 my $smtpaddr = (($hash->{'smtp'} =~ /^(\d{1,3})\.(\d{1,3})\.(\d{1,3})\.(\d{1,3})$/) ?
             pack('C4',$1,$2,$3,$4) :
             (gethostbyname($hash->{'smtp'}))[4]);

 if (!$smtpaddr) { return "Host not found: ".$hash->{'smtp'}; }

 my $s = &FileHandle::new(FileHandle);
 if (!socket($s, AF_INET, SOCK_STREAM, $proto)) {
    return "Socket creation failed!" }

 if (!connect($s, pack('Sna4x8', AF_INET, $port, $smtpaddr))) {
   return "Socket connection failed"; }

 my($oldfh) = select($s); $| = 1; select($oldfh);

 $_ = <$s>; if (/^[45]/) { close $s; return "Server not available "; }
 print $s "helo $client\r\n";
 $_ = <$s>; if (/^[45]/) { close $s; return "Communication error at helo $client"; }
 
 print $s "mail from: <$hash->{'from'}>\r\n";
 $_ = <$s>; if (/^[45]/) { close $s; return "Communication error at 'mail from: $hash->{'from'}'"; }

 { local $^W;
  foreach (split(/, */, $hash->{'to'}),split(/, */, $hash->{'cc'})) {
    next unless /\w/; # a basic sanity check
    if (/<(.*)>/) {
       print $s "rcpt to: $1\r\n";
       } 
      else {
       print $s "rcpt to: <$_>\r\n";
       }
    $_ = <$s>; if (/^[45]/) { close $s; return "User unknown : $hash->{'to'}, [ $hash->{'smtp'} ] "; }
    } #for
 } #block

 print $s "data\r\n";
 $_ = <$s>; if (/^[45]/) { close $s; return "Communication error at command 'data':\n$_"; }

 print $s "To: $hash->{to}\r\n";
 print $s "From: $hash->{from}\r\n";
 print $s "Cc: $hash->{cc}\r\n" if defined $hash->{cc};
 print $s "Reply-to: $hash->{replyto}\r\n" if $hash->{replyto};
 print $s "X-Mailer: Simple perl mailer\r\n"; #this may be omitted
 print $s "Subject: $hash->{subj}\r\n\r\n";
 my $body = $hash->{body};
 $body.=$file; #"attach" file
 $body =~ s/(\A|[^\r])\n/$1\r\n/sg; #ensure each line ends with \r\n
 $body =~ s/\r\n\.\r\n/\r\n\*\r\n/sg;
 print $s $body;
 #end transmission:  
 print $s "\r\n.\r\n";
 $_ = <$s>; if (/^[45]\d* (.*)$/) { close $s; return "Body transmission failed:\n$_\n"; }
 print $s "quit\r\n";
 $_ = <$s>;
 close $s;
 return 0;
}

