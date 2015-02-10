package MALIGN::PAML;

use strict;
use warnings;
use Carp;
use Cwd;

use MALIGN::Alignment;
use MALIGN::Alignment_writer;

my $TMPDIR = $ENV{TMPDIR} || "/tmp";

####
sub run_yn00 {
	my ($alignment_obj) = @_;

	unless (ref $alignment_obj) {
		confess "Error, need alignment_obj for parameter.";
	}
	
	## get the alignment_text
	my $alignment_text = &MALIGN::Alignment_writer::get_alignment_text($alignment_obj, "phylip"); 

	my $workdir = cwd();

	## establish a temp area:
	my $curr_tmpdir = "$TMPDIR/__pamlpm$$." . $ENV{HOSTNAME};
	
	mkdir ($curr_tmpdir) or die "Error, cannot create directory $curr_tmpdir";
	chdir $curr_tmpdir or die "Error, cannot cd to $curr_tmpdir";

	my $alnfile = "$curr_tmpdir/tmp.pamlpm.$$";
	# write the aln file:
	{
		open (my $fh, ">$alnfile") or die "Error, cannot write to $alnfile";
		print $fh $alignment_text;
		close $fh;
	}

	## wrote the yn00.ctl file

	my $yn00ctl_template = <<__EOTEMPLATE;

      seqfile = __INPUTFILE__ * sequence data file name
      outfile = __OUTPUTFILE__   * main result file
      verbose = 0  * 1: detailed output (list sequences), 0: concise output

	  icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below

    weighting = 0  * weighting pathways between codons (0/1)?
   commonf3x4 = 0  * use one set of codon freqs for all pairs (0/1)? 
*       ndata = 1


* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,
* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., 
* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., 
* 10: blepharisma nu.
* These codes correspond to transl_table 1 to 11 of GENEBANK.

__EOTEMPLATE

;
	my $outfile = "yn";
	
	$yn00ctl_template =~ s/__INPUTFILE__/$alnfile/;
	$yn00ctl_template =~ s/__OUTPUTFILE__/$outfile/;

	{ # write the yn00.ctl file
		open (my $fh, ">yn00.ctl") or die "Error, cannot write to yn00.ctl";
		print $fh $yn00ctl_template;
		close $fh;
	}

	my $cmd = "yn00 1>/dev/null 2>/dev/null";
	my $ret = system $cmd;
	if ($ret) {
		die "Error, cmd $cmd died with ret, at $curr_tmpdir";
	}
	
	my @files = (<*>);
	my %outputs;
	foreach my $file (@files) {
		$outputs{$file} = `cat $file`;
	}

	## go back to working directory
	chdir ($workdir) or die "Error, cannot cd back to $workdir";

	unless ($outputs{yn} =~ /\w/) {
		## didn't work, and PAML fails to return nonzero on death for some reason...
		confess "Error, PAML yn00 failed to produce results.  See $curr_tmpdir ";
	}
		
	# cleanup the tmp area:
	`rm -rf $curr_tmpdir`;
	
	if ($?) {
		confess "error, unable to remove the tmp files at $curr_tmpdir";
	}
	
	## parse the resulting {yn} output using PAML::yn00_output_parser

	return (%outputs);
}



1;


	
