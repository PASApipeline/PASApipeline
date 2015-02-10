#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling);
use POSIX qw(ceil floor);

my $usage = <<__EOUSAGE__;

############################################################################
#
#   usage: $0 -A 1,10,2  -B 1,10,2 < val_pairs_file
#
############################################################################
#
#  --pairs <string>
#
#  --col_A_binning || -A <string>        min,max,step
#
#  --col_B_binning || -B <string>        min,max,step
#
#  Optional:
#
#  --log2_data                         bin log2 of value instead of value itself. 
#  --log2_z                            plot the log(bin_val)
#
#  --plot_contour                 writes \$pairs.contour.pdf file
#  --swap_axes
#
############################################################################


Then, plot using R:
 
  (example, eventually do this automatically)

data = read.table("data.txt", header=T, row.names=1)
z = as.matrix(data)
filled.contour(x=seq(0,15,length.out=nrow(z)), y=seq(0,15,length.out=ncol(z)), z=z, color.palette=terrain.colors,asp=1)

############################################################################


__EOUSAGE__

	;


my ($help_flag, $col_A_binning, $col_B_binning, $take_log);
my $pairs_file;
my $plot_contour = 0;
my $take_log_z;
my $swap_axes;

&GetOptions ( 'h' => \$help_flag,
			  
              'pairs=s' => \$pairs_file,
              'col_A_binning|A=s' => \$col_A_binning,
			  'col_B_binning|B=s' => \$col_B_binning,
			  'log2_data' => \$take_log,
              'log2_z' => \$take_log_z,
              'plot_contour' => \$plot_contour,
              'swap_axes' => \$swap_axes,
              );


if ($help_flag) { die $usage; }

unless ($pairs_file) { die $usage; }
unless ($col_A_binning && $col_B_binning) { die $usage; }

if ($swap_axes) {
    ($col_A_binning, $col_B_binning) = ($col_B_binning, $col_A_binning);
}

my ($col_A_min, $col_A_max, $col_A_step) = split (/,/, $col_A_binning);
my ($col_B_min, $col_B_max, $col_B_step) = split (/,/, $col_B_binning);

unless ($col_A_min < $col_A_max && $col_B_min < $col_B_max 
		&& $col_A_step > 0 && $col_B_step > 0) {
	die "Error, params specified incorrectly.  Please check.  \n";
}

my @matrix2D;


open (my $fh, $pairs_file) or die "Error, cannot open file $pairs_file";
while (<$fh>) {
	my $line = $_;
	chomp;
	my ($valA, $valB, @trash) = split (/\s+/);
	unless ($valA =~ /\d/ && $valB =~ /\d/) {
		print STDERR "Error, do not recognize line:\n$line\n";
	}
	
    if ($swap_axes) {
        ($valB, $valA) = ($valA, $valB);
    }


	if ($take_log) {
		$valA = log($valA)/log(2) unless ($valA == 0);
		$valB = log($valB)/log(2) unless ($valB == 0);
	}
	
	## check values for within range.
	my $out_of_range_flag = 0;
	unless (&is_between($valA, $col_A_min-$col_A_step, $col_A_max)) { 
		#print STDERR "$valA of colA is out of range.\n";
		$out_of_range_flag++;
	}
	unless (&is_between($valB, $col_B_min-$col_B_step, $col_B_max)) {
		#print STDERR "$valB of colB is out of range.\n";
		$out_of_range_flag++;
	}
	if ($out_of_range_flag) {
		#print STDERR "Out of range, skipping line:\n$line\n";
		next;
	}
	
	## compute bins.
	my $binA = ceil( ($valA - $col_A_min) / $col_A_step);
	my $binB = ceil( ($valB - $col_B_min) / $col_B_step);

	$matrix2D[$binA][$binB]++;
}

######  report matrix of counts ######
# print header
#print "_x_";
my $matrix_text = "";
for (my $i = $col_A_min; $i <= $col_A_max; $i += $col_A_step) {
	$matrix_text .= "\t$i";
}
$matrix_text .= "\n";

for (my $j = $col_B_min; $j <= $col_B_max; $j += $col_B_step) {
    $matrix_text .= "$j";
	my $bin_j = ceil (($j - $col_B_min) / $col_B_step);
	for (my $i = $col_A_min; $i <= $col_A_max; $i += $col_A_step) {
		my $bin_i = ceil(($i - $col_A_min) / $col_A_step);
		
		my $val = $matrix2D[$bin_i][$bin_j] || 0;
        $matrix_text .= "\t$val";
	}
    $matrix_text .= "\n";
}

if ($plot_contour) {
    my $matrix_file = "$pairs_file.matrix";
    open (my $ofh, ">$matrix_file") or die "Error, cannot write to $matrix_file";
    print $ofh $matrix_text;
    close $ofh;

    &make_contour_plot($matrix_file);
}
else {
    print $matrix_text;
}


exit(0);



####
sub is_between {
	my ($val, $low, $high) = @_;

	if ($val > $low && $val <= $high) { ## bounds tested are '(..]', particularly important for capturing first bin.  
		return(1);
	}
		
	return(0);
}


####
sub make_contour_plot {
    my ($matrix_file) = @_;

    my $R_script = "$matrix_file.R";
    open (my $ofh, ">$R_script") or die "Error, cannot write to file $R_script";
    
    print $ofh "pdf(file=\"$matrix_file.pdf\")\n";
    print $ofh "data = read.table(\"$matrix_file\")\n";
    print $ofh "z = as.matrix(data)\n";
    
    if ($take_log_z) {
        print $ofh "z = log2(z+1)\n";
    }
    
    print $ofh "filled.contour(x=seq($col_B_min,$col_B_max,length.out=nrow(z)), y=seq($col_A_min,$col_A_max,length.out=ncol(z)), z=z, color.palette=terrain.colors, zlim=c(0,quantile(z, c(0.99))))\n";
    print $ofh "dev.off()\n";
    close $ofh;


    my $cmd = "R --vanilla -q < $R_script";
    my $ret = system($cmd);
    
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }
    
    return;
}
    
