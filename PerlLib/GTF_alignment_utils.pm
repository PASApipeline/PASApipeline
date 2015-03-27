#!/usr/bin/env perl

package GTF_alignment_utils;

use strict;
use warnings;
use Carp;

use Gene_obj;
use Gene_obj_indexer;
use CDNA::Alignment_segment;
use CDNA::CDNA_alignment;
use File::Basename;

__run_test() unless caller;
    


sub index_alignment_objs {
    my ($gtf_alignment_file, $genome_alignment_indexer_href) = @_;
    
    unless ($gtf_alignment_file && -s $gtf_alignment_file) {
        confess "Error, cannot find or open file $gtf_alignment_file";
    }
    unless (ref $genome_alignment_indexer_href) {
        confess "Error, need genome indexer href as param ";
    }

    
    	
	my %genome_trans_to_alignment_segments;
	
    my %trans_to_gene_id;

	open (my $fh, $gtf_alignment_file) or die "Error, cannot open file $gtf_alignment_file";
	while (<$fh>) {
		chomp;
		if (/^\#/) { next; }
		unless (/\w/) { next; }
		
		my @x = split(/\t/);

		unless (scalar (@x) >= 8 && $x[8] =~ /transcript_id/) {
			print STDERR "ignoring line: $_\n";
			next;
		}
		
		my $scaff = $x[0];
		my $type = $x[2];
		my $lend = $x[3];
		my $rend = $x[4];
        my $per_id = $x[5];
        
        unless ($type eq 'exon') { next; }
        

        if ($per_id eq ".") { $per_id = 100; } # making an assumption here.
        
		my $orient = $x[6];
		
		my $info = $x[8];
		
		my @parts = split(/;/, $info);
		my %atts;
		foreach my $part (@parts) {
			$part =~ s/^\s+|\s+$//;
			$part =~ s/\"//g;
			my ($att, $val) = split(/\s+/, $part);
			
			if (exists $atts{$att}) {
				die "Error, already defined attribute $att in $_";
			}
			
			$atts{$att} = $val;
		}
        
		my $gene_id = $atts{gene_id} or die "Error, no gene_id at $_";
		my $trans_id = $atts{transcript_id} or die "Error, no trans_id at $_";
        
		        
        $trans_to_gene_id{$trans_id} = $gene_id;

        push (@{$genome_trans_to_alignment_segments{$scaff}->{$trans_id}}, [$lend, $rend, $orient]);
                
    }
    
    
    my %scaff_to_align_list;
    
    
	## Output genes in gtf format:

	foreach my $scaff (sort keys %genome_trans_to_alignment_segments) {
        
        my @alignment_accs = keys %{$genome_trans_to_alignment_segments{$scaff}};

        foreach my $alignment_acc (@alignment_accs) {

            my @segments = @{$genome_trans_to_alignment_segments{$scaff}->{$alignment_acc}};
            @segments = sort {$a->[0]<=>$b->[0]} @segments;
            my $orient = $segments[0]->[2];
            if ($orient eq '-') { 
                @segments = reverse @segments;
            }
            my @cdna_align_segs;
            my @coords;
            my $curr_cdna_len = 0;
            foreach my $segment (@segments) {
                my ($lend, $rend, $orient) = @$segment;
                my ($m_lend, $m_rend) = ($curr_cdna_len + 1, $curr_cdna_len + abs($rend - $lend) + 1);
                $curr_cdna_len = $m_rend;
                if ($orient eq '-') {
                    ($m_lend, $m_rend) = ($m_rend, $m_lend);
                }
                my $alignment_segment = new CDNA::Alignment_segment($lend, $rend, $m_lend, $m_rend, 100);
                push (@cdna_align_segs, $alignment_segment);
                
            }
            
            my $cdna_alignment_obj = new CDNA::CDNA_alignment($curr_cdna_len, \@cdna_align_segs);
            $cdna_alignment_obj->set_acc($alignment_acc);
            $cdna_alignment_obj->{genome_acc} = $scaff;

            $cdna_alignment_obj->{gene_id} = $trans_to_gene_id{$alignment_acc};
            
            my $spliced_orient = $orient;
            if ($spliced_orient !~ /^[\+\-]$/) {
                $spliced_orient = '?';
            }
            $cdna_alignment_obj->set_spliced_orientation($spliced_orient);

            $cdna_alignment_obj->{source} = basename($gtf_alignment_file);
            
            if (ref $genome_alignment_indexer_href eq "Gene_obj_indexer") {
                $genome_alignment_indexer_href->store_gene($alignment_acc, $cdna_alignment_obj);
            }
            else {
                $genome_alignment_indexer_href->{$alignment_acc} = $cdna_alignment_obj;
            }
            
            push (@{$scaff_to_align_list{$scaff}}, $alignment_acc);
        }
    }

    return(%scaff_to_align_list);
}


###########
# Testing
##########

sub __run_test {
    
    my $usage = "usage: $0 file.gtf\n\n";

    my $gtf_file = $ARGV[0] or die $usage;

    my $indexer = {};
    my %scaff_to_alignments = &index_alignment_objs($gtf_file, $indexer);
    
    
    foreach my $scaffold (keys %scaff_to_alignments) {

        my @align_ids = @{$scaff_to_alignments{$scaffold}};

        foreach my $align_id (@align_ids) {
            my $cdna_obj = $indexer->{$align_id};

            print $cdna_obj->toString();
        }
    }

    
    
    exit(0);
    


}





1; #EOM

