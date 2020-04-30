#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

# I use ActiveState Perl, which has a much nice package manager called ppm
# ppm install Bio::DB::Sam
# documentation at https://metacpan.org/pod/Bio::DB::Sam
use Bio::DB::Sam;

# https://github.com/MullinsLab/Bio-Cigar
# perl -MCPAN -e shell
# install Bio::Cigar
use Bio::Cigar;

my $usage = "Usage: $0 <infile.bam> <infile.fa>\n";
my $bam = shift or die $usage;
my $fasta = shift or die $usage;

my $sam = Bio::DB::Sam->new(
   -bam   => $bam,
   -fasta => $fasta
);

my @ref = $sam->seq_ids;
my @alignments = $sam->get_features_by_location(-seq_id => $ref[0],
                                                -start  => 0,
                                                -end    => 20000);

# header information
my $bam_object = Bio::DB::Bam->open($bam);
my $header = $bam_object->header;
my $text = $header->text;

for my $a (@alignments) {

   # alignment line
   my $line = $a->tam_line;
   # print "$line\n";

   # get NM tag
   my $nm = $a->get_tag_values('NM');
 
   # alignment information
   my $read_id =  $a->display_name;
   my $seqid  = $a->seq_id;
   my $start  = $a->start;
   my $end    = $a->end;
   my $strand = $a->strand;
   my $mapping_qual = $a->qual;
   if ($strand == 1){
      $strand = "+";
   } elsif ($strand == -1){
      $strand = "-";
   }
   print join("\t", $read_id, $seqid, $start, $end, $strand, $mapping_qual, $nm), "\n";

   # sequence information
   my $ref_dna   = $a->dna;
   my $query_dna = $a->query->dna;
   my @scores    = $a->qscore;
   # print join("\t", $ref_dna, $query_dna, "@scores"), "\n";

   # query sequence information
   my $query_start = $a->query->start;     
   my $query_end   = $a->query->end;
   # print join("\t", $query_start, $query_end, $query_dna), "\n";

   # CIGAR string
   my $cigar  = $a->cigar_str;
   my $biocigar = Bio::Cigar->new($cigar);
   my $ref_len = $biocigar->reference_length;
   my $query_len = $biocigar->query_length;
   # print join("\t", $query_len, $ref_len), "\n";

   # this part is relevant for spliced sequences
   my $j = 1;
   for (my $i = $start; $i <= $end; ++$i){
      # the rpos_to_qpos() function converts reference coordinates to query coordinates
      my ($qpos, $op) = $biocigar->rpos_to_qpos($j);
      # print "$qpos, $op\n";
      ++$j;
   }

}

exit(0);

