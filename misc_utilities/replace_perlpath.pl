#!/usr/bin/env perl

use strict;
use warnings;

## read list of filenames

my $perl_path = $ARGV[0] || "#!/usr/bin/env perl";



while (my $file = <STDIN>) {
	chomp $file;
	
    if (-f $file) {
        my $tmpfile = "$file.tmp.$$";
        open (my $fh, ">$tmpfile") or die "Error, cannot open file $tmpfile";
        open (my $infh, $file) or die "Error, cannot open file $file";
        my $is_header = 1;
		while (<$infh>) {
            if ($is_header && /^\#/ && /perl/) {  # only update those headers that contain perl in them.
                print $fh "$perl_path\n";
            }
            else {
                print $fh $_;
            }
			$is_header = 0; # not the first line any more.
        }
        close $fh;
        close $infh;

        rename ($tmpfile, $file) or die "Error, cannot rename $tmpfile to $file";
    }
    print STDERR "updated $file\n";
}

exit(0);

