#!/usr/bin/env perl

use strict;
use warnings;
use Pasa_init;
use Pasa_conf;
use DBI;
use CGI;
#use CGI::Carp qw(fatalsToBrowser);
use Data::Dumper;
use Mysql_connect;
use Ath1_cdnas;
use Gene_obj;
use Gene_cdna_image;
use GD::SVG;
use Sequence_feature_drawer;
use TierFeatures;
use CGI::Pretty ":standard";
use Pasa_CGI;
use CDNA::Splice_graph_assembler;
use Splice_graph_illustrator;
use Gene_illustrator;
use SequenceTickerIllustrator;

my $cgi = new CGI();
#print $cgi->header('text/html');

my $db = $cgi->param('db');
my $cdna_acc = $cgi->param('cdna_acc');
my $DEBUG = $cgi->param('DEBUG');

unless ($cdna_acc && $db) {
    die "Need db, cdna_acc\n";
}


my $pid = $$;
unless ($ENV{WEBSERVER_TMP}) {
    $ENV{WEBSERVER_TMP} = "/tmp";
}

# some image settings:
my $image_width = 700;
my $ticker_height = 100;
my $splice_graph_height = 100;
my $height_per_pasa_assembly = 20;
my $spacing_between_graph_and_assemblies = 50;
my $spacing_between_pasa_assemblies = 20;
my $left_margin = 50;
my $right_margin = 200;
my $top_margin = 50;
my $bottom_margin = 50;
my $text_spacing = 20;

my $show_all_flag = $cgi->param('SHOW_ALL');
my $show_alignments_flag = $cgi->param('SHOW_ALIGNMENTS');

my $mysql_server = &Pasa_conf::getParam("MYSQLSERVER");
my $mysql_ro_user = &Pasa_conf::getParam("MYSQL_RO_USER");
my $mysql_ro_password = &Pasa_conf::getParam("MYSQL_RO_PASSWORD");

my ($dbproc) = &connect_to_db($mysql_server,$db,$mysql_ro_user,$mysql_ro_password);

## get the subcluster
my $query = "select subcluster_id from subcluster_link where cdna_acc = \"$cdna_acc\"";
my $subcluster_id = &very_first_result_sql($dbproc, $query);

## get all the pasa assemblies for that subcluster.
$query = "select cdna_acc from subcluster_link where subcluster_id = $subcluster_id";
my @pasa_asmbls;
my @results = &do_sql_2D($dbproc, $query);
foreach my $result (@results) {
    my ($pasa_acc) = @$result;
    push (@pasa_asmbls, $pasa_acc);
}

## get the alignment objects
my @pasa_alignment_objs;
foreach my $pasa_acc (@pasa_asmbls) {
    my $align_id = &Ath1_cdnas::get_validating_align_id_via_acc($dbproc, $pasa_acc);
    my $alignment_obj = &Ath1_cdnas::create_alignment_obj($dbproc, $align_id);
    push (@pasa_alignment_objs, $alignment_obj);
}

my $splice_graph = new CDNA::Splice_graph_assembler();
$splice_graph->build_splicing_graph(@pasa_alignment_objs);


# print $splice_graph->toString();
my @nodes = $splice_graph->get_graph_nodes();

my $num_pasa_assemblies = scalar (@pasa_alignment_objs);

## Draw the splice graph.
my $image_height = $top_margin + $ticker_height + $splice_graph_height + $spacing_between_graph_and_assemblies
    + ($spacing_between_pasa_assemblies * ($num_pasa_assemblies -1) )
    + ($height_per_pasa_assembly * $num_pasa_assemblies)
    + ($bottom_margin);

my $image = new GD::SVG::Image ($image_width, $image_height);

## set up some colors:
my $black = $image->colorAllocate(0,0,0); # black
my $blue = $image->colorAllocate(0, 0, 255); # blue
my $red = $image->colorAllocate(255, 0, 0); # red

my $exon_color = $black;
my $intron_color = $blue;

## set up the coordinate converter.
# get range of coordinates:
my @all_coords;
foreach my $node (@nodes) {
    push (@all_coords, $node->get_coords());
}
@all_coords = sort {$a<=>$b} @all_coords;
my $lend_coord = shift @all_coords;
my $rend_coord = pop @all_coords;
@all_coords = (); # free it.

my $splice_graph_illustrator = new Splice_graph_illustrator();


my ($box_lend, $box_rend) = ($left_margin, $image_width - $right_margin);
my ($box_top, $box_bottom) =  ($top_margin, $top_margin + $splice_graph_height);

my $coord_converter = new Coordinate_draw_converter([$box_lend, $box_rend],
                                                    [$lend_coord, $rend_coord]);



$splice_graph_illustrator->draw_splice_graph($image, \@nodes, $coord_converter, 
                                             
                                             # region of image to draw the splicing graph
                                             {lend => $box_lend,
                                              rend => $box_rend,
                                              top => $box_top,
                                              bottom => $box_bottom },
                                             
                                             # color prefs
                                             { exon => $black,
                                               intron => $blue, },
                                             
                                             );

## Now, draw each pasa-assembly based gene object, with longest ORF allowing 5' and 3' partials.
my $y_top = $box_bottom;
$y_top += $spacing_between_graph_and_assemblies;

foreach my $pasa_acc (@pasa_asmbls) {
    my $gene_obj = &Ath1_cdnas::get_asmbl_gene_obj($dbproc, $pasa_acc, 1, 1);
    
    
    my $y_bottom = $y_top + $height_per_pasa_assembly;
    
    my $gene_illustrator = new Gene_illustrator();
    $gene_illustrator->draw_gene($image, $gene_obj, $coord_converter, 
                                 
                                 # region to draw the gene
                                 { lend => $box_lend,
                                   rend => $box_rend,
                                   top => $y_top,
                                   bottom => $y_bottom, },
                                 
                                 # color preferences
                                 { exon => $black,
                                   CDS => $red, },
                                                                  
                                 );
    
    ## add the accession text
    $image->string(gdSmallFont, $box_rend + $text_spacing, $y_top, $pasa_acc, $black);

    $y_top = $y_bottom + $spacing_between_pasa_assemblies;
}

# add the ticker
my $ticker_illustrator = new SequenceTickerIllustrator( { major_unit_height => 20,
                                                          minor_unit_height => 10 });

## make a new coord converter, start at pos 0.
my $adj_lend = 0;
my $adj_rend = $rend_coord - $lend_coord -1;

$coord_converter = new Coordinate_draw_converter([$box_lend, $box_rend],
                                                 [$adj_lend, $adj_rend]);

$ticker_illustrator->draw($image, $adj_lend, $adj_rend, 
                          # bounding box
                          { lend => $box_lend,
                            rend => $box_rend,
                            top => $y_top,
                            bottom => $y_top + $ticker_height, },

                          $coord_converter,

                          # color prefs
                          { color => $black }
                          
                          );


my $svg_xml = $image->svg();
$svg_xml =~ s/font=/font-face=/g;
print $svg_xml;


exit(0);


