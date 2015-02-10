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
#   --mclclusters  tribemcl output file
#   --annot_inx    gene structures inx file
#
################################################################

_EOH_

    ;

my ($mclclusters, $annot_inx);

&GetOptions ( "mclclusters=s" => \$mclclusters,
              "annot_inx=s" => \$annot_inx,
              );

unless ($mclclusters && $annot_inx) {
    die $usage;
}



my $gene_indexer = new Gene_obj_indexer( { "use" => $annot_inx } );

main: {

    my @fam_cluster;
    my $curr_fam = undef;
    
    open (my $fh, $mclclusters) or die "Error, cannot open file $mclclusters";
    while (<$fh>) {
        unless (/\w/) { next; }
        chomp;
        my ($fam_id, $gene_id, $rest) = split (/\t/);
        
        if (defined($curr_fam) && $fam_id ne $curr_fam) {
            &process_fam_cluster(@fam_cluster);
            @fam_cluster = ();
        }
        push (@fam_cluster, [$fam_id, $gene_id, $rest] );

        $curr_fam = $fam_id;
    }
    close $fh;

    &process_fam_cluster(@fam_cluster); # get last one.

    exit(0);
}

#### 
sub process_fam_cluster {
    my @fam_cluster = @_;

    my $FAM_ID;

    my @structs;
    foreach my $fam_cluster (@fam_cluster) {
        my ($fam_id, $gene_id, $rest) = @$fam_cluster;
        
        $FAM_ID = $fam_id;
        
        eval {
            # get gene coords:
            my $gene_obj = $gene_indexer->get_gene($gene_id);
            
            my ($contig, $strand, $lend, $rend) = ($gene_obj->{asmbl_id}, $gene_obj->{strand}, sort {$a<=>$b} $gene_obj->get_coords());
            
            push (@structs, { gene_id => $gene_id,
                              rest => $rest,
                              contig => $contig,
                              lend => $lend,
                              rend => $rend,
                              strand => $strand,
                          });
        }
        
    };
    
    if ($@) {
        print STDERR "Error: $@";
    }
    
    my $num_entries = scalar (@structs);
    @structs = sort {$a->{contig} cmp $b->{contig} 
                     ||
                         $a->{lend} <=> $b->{lend}} @structs;
    
    

    print "// FAM: $FAM_ID\tSIZE: $num_entries\n";
    foreach my $struct (@structs) {
        my ($contig, $lend, $rend, $strand, $gene_id, $rest) = ($struct->{contig},
                                                                $struct->{lend},
                                                                $struct->{rend},
                                                                $struct->{strand},
                                                                $struct->{gene_id},
                                                                $struct->{rest});

        print "$FAM_ID\t$contig\t$lend\t$rend\t$strand\t$gene_id\t$rest\n";
    }

    print "\n"; # add spacer


    return;
}




