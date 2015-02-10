package Rlite;

## light interface to some useful R functions.

use strict;
use warnings;
use Carp;

use Resource;

#### 
sub G_test_matrix {
	my ($data_rows_aref, $correction) = @_;  #  ( [ [1, 5], [11, 14], ...])  contingency table, one row_aref at a time.
	
	my @data_rows = @$data_rows_aref;

	## correction can be "williams" or "yates"

	if ($correction && ! ($correction eq 'williams' || $correction eq 'yates')) {
		confess "invalid correction specification";
	}
	
	my $Rsub_file = "Rlite/Rsubs/g.test.r";
	my $Rsub_dir = &Resource::find_resource($Rsub_file) or confess "Error, cannot lcoate file: $Rsub_file";

	my $ncol = scalar @{$data_rows[0]};
	
	## build the Rscript
	my $rscript = "source (\"$Rsub_dir/$Rsub_file\")\n"
		. "x = c(";
	my @vals;
	foreach my $row (@data_rows) {
		push (@vals, @$row);
	}
	$rscript .= join (",", @vals) . ")\n";
	
	$rscript .= "m = matrix(x, ncol=$ncol, byrow=TRUE)\n";
	
	$rscript .= "result=g.test(m";
	if ($correction) {
		$rscript .= " ,correct=\"$correction\"";
	}
	$rscript .= ")\n";
	
	$rscript .= "print(result\$statistic)\n";
	$rscript .= "print(result\$p.value)\n";
	
	my $result = &_Run_Rscript($rscript);
	
	my @txt_rows = split (/\n/, $result);

	pop @txt_rows;
	my $pvalue = pop @txt_rows;
	pop @txt_rows;
	pop @txt_rows;
	my $Gval = pop @txt_rows;
	
	$Gval =~ s/^\s+//;

	return ($Gval, $pvalue);
	
}



################################################
## private:

####
sub _Run_Rscript {
	my($rscript) = @_;

	my $tmpfile = "R.tmp.$$.script";
	open (my $fh, ">$tmpfile") or confess "Error, cannot write tmpfile $rscript";
	print $fh $rscript;
	close $fh;

	## run R, catpure the output
	
	my $command = "R --vanilla < $tmpfile";
	my $result = `$command`;

	if ($?) {
		confess "Error running R command:\n$command\nret($?) ";
	}
	
	unlink($tmpfile);

	return ($result);

}





1;

