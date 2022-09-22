#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;


## we delete all files we don't need in this directory. Be careful in case users try running it somewhere else, outside this dir.
chdir $FindBin::Bin or die "error, cannot cd to $FindBin::Bin";



my @files_to_keep = qw (genome_sample.fasta.gz
						run_sample_pipeline.pl
						all_transcripts.fasta
						FL_accs.txt
						orig_annotations_sample.gff3
						alignAssembly.config
						annotCompare.config
						all_transcripts.fasta.clean
						all_transcripts.fasta.cln
						cleanme.pl
						clusters_of_valid_alignments.txt.gz
						pasa_sample_db.mysqldump.gz
                        __runMe_no_TDecoder.pl
                        gene_gff3_to_introns.sh
                        extract_introns_from_pasa_assemblies.sh
						tdn.accs
                        custom_alignments.gff3
orig_annotations_sample.gtf
runMe.SQLite.sh
runMe.MySQL.sh
__run_sample_pipeline.pl
Makefile
runMe.SQLite.only_stringtie_inputs.pl
);

my %keep = map { + $_ => 1 } @files_to_keep;

if (-d "gmap_db_dir") {
	print STDERR "-removing gmap_db_dir\n";
	`rm -rf gmap_db_dir`;
}

`rm -rf assemblies`; 
`rm -rf pasa_run.log.dir`;
`rm -rf transdecoder.tmp.*`;
`rm -rf blat_out_dir/`;
`rm -rf compreh_init_build`;
`rm -rf __tmp_classify_alt_isoforms`;
`rm -rf __pasa_sample_mydb_pasa.sqlite_SQLite_chkpts`;
`rm -rf sample_mydb_pasa.sqlite.assemblies.fasta.transdecoder_dir.__checkpoints`;
`rm -rf sample_mydb_pasa.sqlite.assemblies.fasta.transdecoder_dir`;
`rm -rf genome_sample.fasta.gmap`;
`rm -rf sample_mydb_pasa.assemblies.fasta.transdecoder_dir`;
`rm -rf sample_mydb_pasa.assemblies.fasta.transdecoder_dir.__checkpoints`;
`rm -rf __pasa_sample_mydb_pasa_mysql_chkpts`;
`rm -rf __chkpts_sample_mydb_pasa.sqlite`;
`rm -rf ./all_transcripts.fasta.clean.transdecoder_dir*`;
`rm -rf ./__chkpts_sample_mydb_pasa`;
`rm -rf ./sample_mydb_pasa.sqlite.assemblies.fasta.transdecoder_dir.__checkpoints_longorfs`;
`rm -rf ./*__checkpoints_longorfs`;
`rm -rf ./pblat_outdir`;
`rm -rf ./genome_sample.fasta.mm2`;

foreach my $file (<*>) {
	
	if (-f $file && ! $keep{$file}) {
		print STDERR "-removing file: $file\n";
		unlink($file);
	}
}


exit(0);
