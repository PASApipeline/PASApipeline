#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;


my %data;

while (<>) {
    if (/^>(\S+)/) {
        if (%data) {
			&report_bed(%data);
		}
		
		%data = ();
		
		my $transcript_acc = $1;
    
		%data = (acc => $transcript_acc,
				 segments => {},  # end5 => end3
			);
		
	}
    elsif (/\s*([\+\-])([^\:]+):(\d+)-(\d+)\s+\((\d+)-(\d+)\)\s+(\d[^\%]+\%)/) {
        #print " $transcript_acc\t$_";
        my $orient = $1;
        my $genome_acc = $2;
        my $genome_end5 = $3;
        my $genome_end3 = $4;
        my $transcript_end5 = $5;
        my $transcript_end3 = $6;
        my $percent_identity = $7;
        
		
		$data{segments}->{$genome_end5} = $genome_end3;
		$data{scaff} = $genome_acc;
	}
}

if (%data) {
	&report_bed(%data);
}

exit(0);


####
sub report_bed {
	my %data = @_;

	my $gene_obj = new Gene_obj();
	
	my $coords_href = $data{segments};
	my $acc = $data{acc};
	my $genome_scaff = $data{scaff};

	if (%$coords_href) {
		
		$gene_obj->populate_gene_object($coords_href, $coords_href);
		
		$gene_obj->{com_name} = $acc;
		$gene_obj->{asmbl_id} = $genome_scaff;
		
		print $gene_obj->to_BED_format();
	}
	
	return;
}
