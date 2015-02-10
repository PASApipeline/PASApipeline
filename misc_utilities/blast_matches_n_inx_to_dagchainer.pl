#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Gene_obj_indexer;
use Carp;    
use Getopt::Long qw(:config no_ignore_case bundling);   

my $usage = <<_EOH_;

################################################################
#
#   --match_file  btab or m2fmt (wu-blast mformat=2  file)
#   --format      btab | m2fmt
#   --query_inx   query inx filel (see index_gff3_file.pl)
#   --search_inx  search database inx file 
#
################################################################

_EOH_

    ;

my ($match_file, $query_inx, $search_inx, $format);

&GetOptions ( "match_file=s" => \$match_file,
              "format=s" => \$format,
			  "query_inx=s" => \$query_inx,
              "search_inx=s" => \$search_inx,
              );

unless ($match_file && $query_inx && $search_inx) {
    die $usage;
}

unless ($format && $format =~ /^(btab|m2fmt)$/) {
	die $usage;
}

main: {
    my $query_gene_indexer = new Gene_obj_indexer( { "use" => $query_inx } );
    my $search_gene_indexer = new Gene_obj_indexer( { "use" => $search_inx } );
    
    open (my $fh, $match_file) or die $!;
    while (<$fh>) {
        chomp;
        my @x = split (/\t/);
        my ($query_acc, $search_acc, $p_value);
		if ($format eq 'btab') {
			($query_acc, $search_acc, $p_value) = ($x[0], $x[5], $x[19]);
		}
		elsif ($format eq 'm2fmt') {
			($query_acc, $search_acc, $p_value) = (@x[0..2]);
		}
		else {
			die "Error, don't understand format: $format";
			# should never get here anyway due to format check at top
		}
		
		eval {
			my $query_gene = $query_gene_indexer->get_gene($query_acc);
			my $search_gene = $search_gene_indexer->get_gene($search_acc);
			
			
			my ($query_mol, $query_end5, $query_end3) = ($query_gene->{asmbl_id}, $query_gene->get_coords());
			my ($search_mol, $search_end5, $search_end3) = ($search_gene->{asmbl_id}, $search_gene->get_coords());
			
			print "$query_mol\t$query_mol:$query_acc\t$query_end5\t$query_end3\t"
				. "$search_mol\t$search_mol:$search_acc\t$search_end5\t$search_end3\t"
				. "$p_value\n";
		};

		if ($@) {
			print STDERR $@;
		}
	}

    exit(0);
}

        





