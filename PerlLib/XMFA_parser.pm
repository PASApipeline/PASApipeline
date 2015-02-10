package XMFA_parser;

use strict;
use warnings;
use File::Basename;
use Carp;
use Data::Dumper;
use Fasta_reader;

### NOTE: This module puts the entire multiple alignment into memory!!! So Beware!  #####


####
sub new {
    my $packagename = shift;
    my ($xmfa_file) = @_;
    
    unless ($xmfa_file && -s $xmfa_file) {
        confess "Error, need xmfa-formatted (mauve-generated) file as parameter";
    }
    
    my $self = { 
        xmfa_file => $xmfa_file,
        
        alignment_blocks => [], # list of alignment block objects
        
        
    };

    bless ($self, $packagename);

    $self->_parse_xmfa();
    
    return($self);
}


####
sub _parse_xmfa {
    my $self = shift;

    my $xmfa_file = $self->{xmfa_file};
    
    
    my %genomes;
        
    my %seqs;
        
    my $done_parsing_header = 0;
    
    my $curr_acc = "";
    
    my %genome_to_scaffold_name_lookup;

    open (my $fh, $xmfa_file) or die "Error, cannot open file $xmfa_file";
    while (<$fh>) {
        #print;
        chomp;
        if (! $done_parsing_header && /^\#/) {
            my @x = split(/\s+/);
            if ($x[0] =~ /Sequence\d+File/) {
                $genomes{$x[1]} = $x[1];
            }
        }
        elsif (/^\>/) {
            my ($carat, $coords, $orient, $genome) = split(/\s+/);
            my $acc = join("$;", $genome, $orient, $coords);
            $curr_acc = $acc;
            
            if (! $done_parsing_header) {
                $done_parsing_header = 1;
                
                $self->_parse_scaffold_names_from_genome_files(\%genome_to_scaffold_name_lookup, \%genomes, dirname($xmfa_file));
                
            }
            
            
        }
        elsif (/^=/) {
            $self->_process_alignment_block(\%seqs, \%genome_to_scaffold_name_lookup);
            # reinit
            %seqs = ();
            $curr_acc = "";
        }
        else {
            $seqs{$curr_acc} .= $_;
        }
    }

    close $fh;

    if (%seqs) {
        # get last one
        $self->_process_alignment_block(\%seqs, \%genome_to_scaffold_name_lookup);
    }
    
    
    
    return;
}


####
sub _process_alignment_block {
    my $self = shift;
    my ($seqs_href, $scaffold_lookup_href) = @_;
    
    
    my @alignments;

    foreach my $acc (keys %$seqs_href) {
        my ($genome, $orient, $coords) = split(/$;/, $acc);
        
        my ($mol, $lend, $rend) = split(/[:\-]/, $coords);
        
        if ($lend == 0 && $rend == 0) { 
            ## no alignment, just all gaps
            ## ignore it
            next;
        }
        
        my $aligned_seq = $seqs_href->{$acc};
        my ($scaffold_name, $scaff_lend, $scaff_rend)  = $self->_concat_scaff_to_indiv_scaff($scaffold_lookup_href->{$genome}, $lend, $rend);
        unless ($scaffold_name) {
            print STDERR Dumper($scaffold_lookup_href);
            confess "Error, cannot find scaffold name for $genome [$mol-1] ";
        }
        
        #print "-adding: $genome, $scaffold_name, $scaff_lend-$scaff_rend, $orient\n";
        
        my $aligned_seq_obj = Aligned_sequence->new($genome, $scaffold_name, $scaff_lend, $scaff_rend, $orient, $aligned_seq);
        
        push (@alignments, $aligned_seq_obj);
        
        
    }
    
    
    if (scalar @alignments >= 2) {
        my $alignment_block = Alignment_block->new(@alignments);
            
        $self->_add_alignment_block($alignment_block);
    }
    
        
    return;
}


####
sub _add_alignment_block {
    my $self = shift;
    my ($alignment_block) = @_;

    unless ($alignment_block) {
        confess "need alignment block as param";
    }

    push (@{$self->{alignment_blocks}}, $alignment_block);

    return;
}


####
sub get_alignment_blocks {
    my $self = shift;
    return(@{$self->{alignment_blocks}});
}

####
sub find_block_containing_seqrange {
    my $self = shift;
    my ($genome, $contig, $lend, $rend) = @_;

    unless (defined($genome) && defined($contig) && defined($lend) && defined($rend)) {
        confess "Error, see required params";
    }

    ($lend, $rend) = sort {$a<=>$b} ($lend, $rend);


    #print "Searching blocks for: [$genome] [$contig] [$lend-$rend]\n";
    
    my $found_alignment_block = undef;

    eval {
        foreach my $alignment_block ($self->get_alignment_blocks()) {
            
            if ($alignment_block->contains_seqrange($genome, $contig, $lend, $rend)) {
                $found_alignment_block = $alignment_block;
                #print "Got it.\n";
                last;
            }
        }
    };

    if ($found_alignment_block) {
        return $found_alignment_block;
    }
    else {
        #print "Block not found.\n";
        return undef;
    }
}


####
sub print() {
    my $self = shift;

    my @alignment_blocks = $self->get_alignment_blocks();
    foreach my $alignment_block (@alignment_blocks) {

        print $alignment_block->toString();
    }
    
    return;
}


####
sub _parse_scaffold_names_from_genome_files {
    my $self = shift;
    
    my ($genome_to_scaff_lookup_href, $genome_to_file_href, $dirname) = @_;
    
    foreach my $genome (keys %$genome_to_file_href) {
        
        my $genome_fasta_file = "$dirname/" . $genome_to_file_href->{$genome};
        print STDERR "// parsing scaffold names from $genome_fasta_file\n";

        my $fasta_reader = new Fasta_reader($genome_fasta_file);

        my $scaffold_lend = 0;

        while (my $seq_obj = $fasta_reader->next()) {

            my $acc = $seq_obj->get_accession();
            my $sequence = $seq_obj->get_sequence();
            
            my $seq_len = length($sequence);
            
            my $scaff_struct = { acc => $acc,
                                 lend => 1,
                                 rend => $seq_len,
                                 scaffold_lend => $scaffold_lend + 1,
                                 scaffold_rend => $scaffold_lend + $seq_len,
                             };

            $scaffold_lend += $seq_len;
            
            push (@{$genome_to_scaff_lookup_href->{$genome}}, $scaff_struct);
            
        }
        
    }
    
    return;
}


sub _concat_scaff_to_indiv_scaff {
    my $self = shift;
    my ($scaff_struct_list_aref, $lend, $rend) = @_;

    foreach my $scaff_struct (@$scaff_struct_list_aref) {
        my ($scaff_name, $scaff_lend, $scaff_rend) = ($scaff_struct->{acc}, $scaff_struct->{scaffold_lend}, $scaff_struct->{scaffold_rend});
        
        if ($lend >= $scaff_lend && $lend <= $scaff_rend) {  ## NOTE: doesn't take into account alignments that may traverse neighboring concatenated contigs

            ## found scaffold

            my $rel_lend = $lend - $scaff_lend + 1;
            my $rel_rend = $rend - $scaff_lend + 1;
            
            return($scaff_name, $rel_lend, $rel_rend);
        }
    }

    print STDERR Dumper($scaff_struct_list_aref);
    
    confess "Error, cannot map region ($lend-$rend) to a scaffold entry.";
}




########################
package Alignment_block;

use strict;
use warnings;

use Carp;
use Storable qw (dclone);


sub new {
    my $packagename = shift;
    my @aligned_seqs = @_;

    unless (scalar @aligned_seqs >= 2) {
        die "Error, need at least two aligned sequence objects as parameters";
    }

    my $self = {
       aligned_seqs => [@aligned_seqs],
    };

    bless($self, $packagename);
    
    return($self);
}

####
sub get_aligned_seqs {
    my $self = shift;

    return(@{$self->{aligned_seqs}});
}


####
sub contains_seqrange {
    my $self = shift;
    my ($genome, $mol, $lend, $rend) = @_;

    unless (defined($genome) && defined($mol) && defined($lend) && defined($rend) ) {
        confess "error, see required params";
    }
    
    
    #print "// block\n";
    foreach my $aligned_seq_obj ($self->get_aligned_seqs()) {

        #if ($aligned_seq_obj->{genome} eq $genome && $aligned_seq_obj->{scaffold} eq $mol) {
        #    print $aligned_seq_obj->header();
        #    if ($aligned_seq_obj->{lend} > $rend) {
        #        print "Not found: $genome, $mol, $lend-$rend\n";
        #        die; # throw exception to be caught
        #    }
        #}
        
        if ($aligned_seq_obj->contains_seqrange($genome, $mol, $lend, $rend)) {
            #print "block:  found seq containing range.\n";
            return(1);
        }
    }

    
    return(0); # doesn't contain region
}

####
sub contains_genome {
    my $self = shift;
    my ($genome_name) = @_;

    unless (defined $genome_name) {
        confess "Error, need genome name as parameter";
    }

    
    foreach my $aligned_seq_obj ($self->get_aligned_seqs()) {

        if ($aligned_seq_obj->{genome} eq $genome_name) {
            return(1);
        }
    }
    
    return(0); # doesn't contain genome
}


####
sub toString {
    my $self = shift;
    
    my @aligned_seqs = $self->get_aligned_seqs();

    my $text = "// Alignment block:\n";

    foreach my $aligned_seq (@aligned_seqs) {

        $text .= $aligned_seq->toString();
    }

    return($text);
}

    
####
sub trim_block_to_genome_coords {
    my $self = shift;
    my ($genome, $scaffold, $lend, $rend) = @_;
    
    my $ref_aligned_seq;
    my @other_aligned_seqs;
    
    foreach my $aligned_seq_obj ($self->get_aligned_seqs()) {
        
        if ($aligned_seq_obj->contains_seqrange($genome, $scaffold, $lend, $rend)) {
            
            if ($ref_aligned_seq) {
                confess "Error, already set ref_aligned_seq, shouldn't encounter this region twice in the same block";
            }
            $ref_aligned_seq = $aligned_seq_obj;
        }
        else {
            
            push (@other_aligned_seqs, $aligned_seq_obj);
        }
        
    }
    
    ## determine the amount of trimming to do:
    my ($left_trim, $right_trim) = $self->_compute_trim_lengths($ref_aligned_seq, $lend, $rend);

    my @trimmed_seqs;
    foreach my $aligned_seq ($ref_aligned_seq, @other_aligned_seqs) {

        my $copy_aligned_seq = dclone($aligned_seq);
        $copy_aligned_seq->trim_alignment($left_trim, $right_trim);
        push (@trimmed_seqs, $copy_aligned_seq);
    }

    my $trimmed_block = Alignment_block->new(@trimmed_seqs);
    
    return($trimmed_block);
        
}


####
sub _compute_trim_lengths {
    my $self = shift;
    my ($aligned_seq, $lend, $rend) = @_;
    
    ($lend, $rend) = sort {$a<=>$b} ($lend, $rend);
    
    my $align_orient = $aligned_seq->{orient};
    my $align_lend = $aligned_seq->{lend};
    my $align_rend = $aligned_seq->{rend};


    #print "Computing trim lengths: starting($align_lend-$align_rend), target: ($lend-$rend)\n";
    
    
    my $alignment_chars = $aligned_seq->{aligned_seq};
    my @chars = split(//, $alignment_chars);
    
    my $left_trim = 0;
    my $right_trim = 0;
    
    if ($align_orient eq '+') {
        my $chars_encountered = 0;
        
        # trim left
        for (my $i = 0; $i <= $#chars; $i++) {
            my $char = $chars[$i];
            if ($char =~ /\w/) {
                $chars_encountered++;
                if ($align_lend + $chars_encountered - 1 == $lend) {
                    $left_trim = $i;
                    last;
                }
            }
        }
        
        # trim right
        $chars_encountered = 0; ## reinit
        for (my $i = $#chars; $i >= 0; $i--) {
            my $char = $chars[$i];
            if ($char =~ /\w/) {
                $chars_encountered++;
                if ($align_rend - $chars_encountered + 1 == $rend) {
                    $right_trim = $#chars - $i;
                    last;
                }
            }
        }

        return($left_trim, $right_trim);
    }
    else {
        ## minus strand
        
        my $chars_encountered = 0;
        
        # trim left
        for (my $i = 0; $i <= $#chars; $i++) {
            my $char = $chars[$i];
            if ($char =~ /\w/) {
                $chars_encountered++;
                if ($align_rend - $chars_encountered + 1 == $rend) {
                    $left_trim = $i;
                    last;
                }
            }
        }
        
        # trim right
        $chars_encountered = 0; ## reinit
        for (my $i = $#chars; $i >= 0; $i--) {
            my $char = $chars[$i];
            if ($char =~ /\w/) {
                $chars_encountered++;
                if ($align_lend + $chars_encountered - 1 == $lend) {
                    $right_trim = $#chars - $i;
                    last;
                }
            }
        }
        
        return($left_trim, $right_trim);
    }
}



####
sub to_pretty_malign_text {
    my $self = shift;
    
    my @aligned_seqs = $self->get_aligned_seqs();
    
    my $REGION_LEN = 60;
    my $curr_pos = 0;
    
    my $align_text = "";
    my $done = 0;
    while (! $done) {


        my @aligns_for_consensus_check;


        foreach my $aligned_seq (@aligned_seqs) {
            
            my $genome = $aligned_seq->{genome};
            my $scaffold = $aligned_seq->{scaffold};
            my $seq_string = $aligned_seq->{aligned_seq};
            
            if (length($seq_string) <= $curr_pos + $REGION_LEN) {
                $done = 1;
            }

            my $acc_name = substr("$genome;$scaffold", 0, 30);
            my $align_region_string = uc substr($seq_string, $curr_pos, $REGION_LEN);

            $align_text .= "$acc_name\t$align_region_string\n";
            
            my @chars = split(//, $align_region_string);
            push (@aligns_for_consensus_check, [@chars]);
        }
        $curr_pos += $REGION_LEN;
        
        my $consensus_text = &build_consensus_line(@aligns_for_consensus_check);
        $align_text .= (" " x 30) . "\t$consensus_text\n";
        $align_text .= "\n"; # spacer
        
        

    }

    return($align_text);
}


sub build_consensus_line {
    my (@aligns) = @_;

    my $consensus_text = "";
    for (my $i = 0; $i <= $#{$aligns[0]}; $i++) {

        my %char_counter;
        foreach my $align (@aligns) {
            $char_counter{ uc $align->[$i] }++;
        }
        my @diff_chars = keys %char_counter;
        if (scalar @diff_chars == 1 && $diff_chars[0] =~ /\w/) {
            $consensus_text .= "*";
        }
        else {
            $consensus_text .= " ";
        }
    }

    return($consensus_text);
}



#########################
package Aligned_sequence;

use strict;
use warnings;

use Carp;


####
sub new {
    my $packagename = shift;
    my ($genome, $scaffold, $lend, $rend, $orient, $aligned_seq) = @_;
    
    unless ($genome && $scaffold && $lend && $rend && $orient && $aligned_seq) {
        confess "Error, see required params";
    }

    my $self = {
        genome => $genome,
        scaffold => $scaffold,
        lend => $lend,
        rend => $rend,
        orient => $orient,
        aligned_seq => $aligned_seq,
    };

    bless ($self, $packagename);

    return($self);
}

####
sub header {
    my $self = shift;

    my $text = join("\t", "#SEQ:", $self->{genome}, $self->{scaffold}, $self->{lend}, $self->{rend}, $self->{orient}) . "\n";
   
    
    return($text);
}


####
sub toString {
    my $self = shift;

    my $text = $self->header() . $self->{aligned_seq} . "\n";
    
    return($text);
}

####
sub contains_seqrange {
    my $self = shift;
    my ($genome, $scaffold, $lend, $rend) = @_; 

    ($lend, $rend) = sort {$a<=>$b} ($lend, $rend); # ensure properly sorted

    #if ($genome eq $self->{genome} && $scaffold eq $self->{scaffold}) {
    #    print "\tchecking: [$scaffold] vs. [$self->{scaffold}], $self->{lend} <= $lend && $rend <= $self->{rend}\n";
    #}
    
    if ($self->{genome} eq $genome
        &&
        $self->{scaffold} eq $scaffold
        &&
        $self->{lend} <= $lend && $rend <= $self->{rend}) {
        
        #print "FOUND IT\n";

        return(1);
    }
    
    #print "not found\n";
    
    return(0); # not within range
}

####
sub trim_alignment {
    my $self = shift;
    my ($left_trim, $right_trim) = @_;

    my $aligned_seq = $self->{aligned_seq};
    
    my $left_over = length($aligned_seq) - $left_trim - $right_trim;
    
    #print "TRIMMING ALIGNMENT($left_trim, $right_trim) from alignment of length: " . length($aligned_seq) . ", leaving: $left_over\n";
    
    my @chars = split(//, $aligned_seq);
    
    my $left_chars_trimmed = 0;
    for (my $i = 0; $i < $left_trim; $i++) {
        my $char = $chars[0];
        if (defined($char) && $char =~ /\w/) {
            $left_chars_trimmed++;
        }
        shift @chars;
    }
    
    my $right_chars_trimmed = 0;
    for (my $i = 0; $i < $right_trim; $i++) {
        my $char = $chars[$#chars];
        if (defined($char) && $char =~ /\w/) {
            $right_chars_trimmed++;
        }
        pop @chars;
    }

    ## update alignment info
    $self->{aligned_seq} = join("", @chars);
    if ($self->{orient} eq '+') {
        $self->{lend} += $left_chars_trimmed;
        $self->{rend} -= $right_chars_trimmed;
    }
    else {
        ## minus orient
        $self->{rend} -= $left_chars_trimmed;
        $self->{lend} += $right_chars_trimmed;
    }


    return;
}



1; #EOM
