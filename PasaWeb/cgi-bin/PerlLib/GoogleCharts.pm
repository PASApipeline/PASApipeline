package GoogleCharts;

use strict;
use warnings;
use Carp;


#########################
##  Barchart ############
#########################



#  Define the column format
#
#   my $columns_aref = [ 
#                         ['string', 'genome'],       # category label (y-axis)
#                         ['number', 'gene count'],   # measurements (x-axis)
#                        ];
# 
#   my $data_aref = [
#                        ['genomeA_name', 500],    # name of genome, count of genes
#                        ['genomeB_name', 250], 
#                   ];
#
#   my $selection_handler_js = # "alert(orgname + \" \" + category);\n" .
#        "window.location.href=\"?project=$project&orthoClusterScan=1&orgname=\" + orgname + \"&category=\" + category;\n";



####
sub draw_barchart {
    my %chart_params = @_;
    
        
    my $columns_aref = $chart_params{columns} or confess "Error, no columns set";
    my $data_aref = $chart_params{data} or confess "Error no data set";
    my $stacked_flag = $chart_params{stacked} || 0;
    my $div_id = $chart_params{div_id} || confess "Error, need div_id set.";
    my $title = $chart_params{title} || $div_id;
    my $selection_handler_js = $chart_params{selection_handler_js} || "";
    my $chart_height = $chart_params{height};
    my $chart_width = $chart_params{width} || 400;
    
    my $barchart_function = "barchart_$div_id";

       
    my $text = "";

    $text .= "<script type=\"text/javascript\" src=\"https://www.google.com/jsapi\"></script>\n";
    $text .= "<script type=\"text/javascript\">\n";
    $text .= "google.load(\"visualization\", \"1\", {packages:[\"corechart\"]});\n";
    $text .= "google.setOnLoadCallback($barchart_function);\n";
    
    $text .= "function $barchart_function\() {\n";
    $text .= "    var data = new google.visualization.DataTable();\n";
    
    my $top_col_name = $columns_aref->[0]->[1];
    
    foreach my $column (@$columns_aref) {
        my ($data_type, $column_name) = @$column;
        $text .= "data.addColumn(\'$data_type\', \'$column_name\');\n";
    }
    
    $text .= "    data.addRows([\n";
    for (my $i = 0; $i <= $#$data_aref; $i++) {
        my ($name, @vals) = @{$data_aref->[$i]};
        
        $text .= "        [\'$name\', " . join(", ", @vals) . "]";
        
        unless ($i == $#$data_aref) {
            $text .= ",\n";
        }
        $text .= "\n";
    }
    $text .= "       ]);\n";
    
    unless ($chart_height) {
        $chart_height = 30 * scalar(@$data_aref);
    }
    
    
    $text .= "    var options = {\n";
    $text .= "        width: $chart_width, height: $chart_height,\n";
    $text .= "        title: \'$title\',\n";
    $text .= "        vAxis: {title: \'$top_col_name\',  titleTextStyle: {color: 'red'}},\n";
    
    if ($stacked_flag) {
        $text .= "        isStacked: true,\n";
    }
    
    $text .= "    };\n";
    
    
    $text .= "var chart = new google.visualization.BarChart(document.getElementById(\'$div_id\'));\n";
    
    $text .= "function processChartSelection_$div_id() {\n"
        . "    s = chart.getSelection();\n"
        . "    r = s[0].row;\n"
        . "    c = s[0].column;\n"
        . "    val = data.getValue(r,c);\n"
        . "    orgname = data.getValue(r,0);\n"
        . "    category = data.getColumnLabel(c);\n"
        #. "    alert(\"row: +\" + r + \", col: +\" + c + \", val: \" + val + \", org: \" + orgname + \", cat: \" + category);\n"
        . "    selection_handler_$div_id(orgname, category);\n"    
        . "}\n";
    
    
    $text .= "google.visualization.events.addListener(chart, 'select', processChartSelection_$div_id);\n";
    
    
    $text .= "function selection_handler_$div_id(orgname, category) {\n"
        . "    $selection_handler_js\n"
        . "}\n\n";
    
    $text .= "    chart.draw(data, options);\n";
    $text .= "}\n";
    $text .= "</script>\n";
    
    $text .= "\n\n<div id=\"$div_id\" style=\"float:left\" ></div>\n\n";
    
    return ($text);
}


##################
## Column Chart ##
##################

=example_column_chart

my %chart_params = ( 

                     column_category_label => 'Year',
                     ordered_column_labels => [ '2004', '2005', '2006', '2007' ],
                     ordered_variable_names => [ 'Sales', 'Expenses' ],
                     data_values => { '2004' => { 'Sales' => 1000,
                                                  'Expenses' => 400,
                                              },
                                      '2005' => { 'Sales' => 1170,
                                                  'Expenses' => 460,
                                              },
                                      '2006' => { 'Sales' => 660,
                                                  'Expenses' => 1120,
                                              },
                                      '2007' => { 'Sales' => 1030,
                                                  'Expenses' => 540,
                                              }
                                  },
                     div_id => "my_div",
                     );

my $chart_text = &draw_columnchart(%chart_params);

=cut


sub draw_columnchart {
    my %chart_params = @_;

    my $column_category_label = $chart_params{column_category_label} or die "Error, need column_category_label";
    my $column_labels_aref = $chart_params{ordered_column_labels} or die "Error, need column_labels";

    my $variables_aref = $chart_params{ordered_variable_names} or die "Error, need ordered_variable_names"; 
    my $data_values_href = $chart_params{data_values} or die "Error, need data values";
    my $div_id = $chart_params{div_id} or die "Error, need div_id";
    my $width = $chart_params{width} || 500;
    my $height = $chart_params{height} || 500;
    my $chart_title = $chart_params{title} || "Chart Title";
    my $y_axis_title = $chart_params{y_axis_title} || "y axis title";  
    

    my $text = "";

    $text .= "<script type=\"text/javascript\" src=\"https://www.google.com/jsapi\"></script>\n";
    $text .= "<script type=\"text/javascript\">\n";
    $text .= "google.load(\"visualization\", \"1\", {packages:[\"corechart\"]});\n";
    $text .= "google.setOnLoadCallback(drawChart_$div_id);\n";
    $text .= "function drawChart_${div_id}\() {\n";
    $text .= "     var data = google.visualization.arrayToDataTable([\n";
    $text .= "        ['" . join("', '", $column_category_label, @$variables_aref) . "'],\n";

    for (my $i = 0; $i <= $#$column_labels_aref; $i++) {
        my $column_label = $column_labels_aref->[$i];
        $text .= "[\'$column_label\'";
        foreach my $var (@$variables_aref) {
            my $val = $data_values_href->{$column_label}->{$var};
            unless (defined $val) {
                die "Error, need value for $column_label -> $var ";
            }
            $text .= ", " . $val;
        }
        $text .= "]";
        if ($i != $#$column_labels_aref) {
            $text .= ",";
        }
        
        $text .= "\n";
    }
    
    # example from google.
    # ['Year', 'Sales', 'Expenses'],
    # ['2004',  1000,      400],
    # ['2005',  1170,      460],
    # ['2006',  660,       1120],
    # ['2007',  1030,      540]
    
    
    $text .= "]);\n";

    $text .= "  var options = {\n";
    $text .= "  title: \'${chart_title}\',\n";
    $text .= "  hAxis: {title: \'${y_axis_title}\', titleTextStyle: {color: 'red'}}\n";
    $text .= "};\n";
    
    $text .= "  var chart = new google.visualization.ColumnChart(document.getElementById('chart_div_$div_id'));\n";
    $text .= "  chart.draw(data, options);\n";
    $text .= " }\n";
    $text .= "</script>\n";
    $text .= "<div id=\"chart_div_$div_id\" style=\"width: ${width}px; height: ${height}px;\"></div>\n";

    return($text);
}



sub draw_piechart {
    my %chart_params = @_;

    my $category_name_description = $chart_params{category_name_description} || "category name description";
    my $category_value_description = $chart_params{category_value_description} || "category value description";
    
    my $data_href = $chart_params{data} or die "Error, need data";
    my $chart_title = $chart_params{title};;
    my $width = $chart_params{width} || 500;
    my $height = $chart_params{height} || 500;
    my $div_id = $chart_params{div_id} or die "Error, need div_id";
    

    my $text = "";
    
    $text .= "<script type=\"text/javascript\" src=\"https://www.google.com/jsapi\"></script>\n";
    $text .= "<script type=\"text/javascript\">\n";
    $text .= "google.load(\"visualization\", \"1\", {packages:[\"corechart\"]});\n";
    $text .= "google.setOnLoadCallback(drawChart);\n";
    $text .= "function drawChart() {\n";
    $text .= "   var data = google.visualization.arrayToDataTable([\n";

        
    # example from google:
    # ['Task', 'Hours per Day'],
    # ['Work',     11],
    # ['Eat',      2],
    # ['Commute',  2],
    # ['Watch TV', 2],
    # ['Sleep',    7]

    $text .= "[ \'$category_name_description\', \'$category_value_description\' ],\n";
    my @data_keys = keys %$data_href;
    for (my $i = 0; $i <= $#data_keys; $i++) {
        my $key = $data_keys[$i];
        my $val = $data_href->{$key};
        $text .= "[ \'$key\', $val ]";
        if ($i != $#data_keys) {
            $text .= ",";
        }
        $text .= "\n";
    }
    
    $text .= "]);\n";

    $text .= "   var options = {\n";
    $text .= "   title: \'$chart_title\'\n";
    $text .= " };\n";

    $text .= "var chart = new google.visualization.PieChart(document.getElementById(\'$div_id\'));\n";
    $text .= "chart.draw(data, options);\n";
    $text .= "}\n";
    $text .= "</script>\n";
  
    $text .= "<div id=\"$div_id\" style=\"width: ${width}px; height: ${height}px;\"></div>\n";
    
    return($text);
}



1; #EOM
