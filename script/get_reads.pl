#!/usr/bin/env perl

use strict;
use warnings;
use Bio::DB::Sam;
use Getopt::Std;
use Parallel::ForkManager;

my %opts = ();
getopts('b:t:l:i:s:p:h:f:', \%opts);

if ( $opts{'h'} ||
   !exists $opts{'b'} ||
   !exists $opts{'l'} ||
   !exists $opts{'f'} ||
   !exists $opts{'i'}
){
   usage();
}

# arguments are explained in the usage at the end of the script
my $bam_file   = $opts{'b'};
my $list       = $opts{'l'};
my $fasta_file = $opts{'f'};
my $idx_file   = $opts{'i'};

my $num_split  = 40;
my $processes  = 8;
my $tmp_dir    = "/tmp/";

if ($opts{'t'}){
   if (!-d $tmp_dir){
      die("$tmp_dir does not exist\n");
   } else {
      warn("Using $tmp_dir as temporary directory\n");
      $tmp_dir = $opts{'t'};
   }
}

if ($opts{'p'}){
   $processes = $opts{'p'};
}

if ($opts{'s'}){
   $num_split = $opts{'s'};
}

# store list of read names
my %list = store_read($list);
my $n = scalar(keys %list);
warn("Stored $n reads\n");

# get chromosome sizes for splitting
my %ref_size = get_ref_size($idx_file);
my @ref = keys %ref_size;

warn("Using $processes processors\n");
my $manager = new Parallel::ForkManager($processes);

# store list of files to merge back together
my @chunk = ();

foreach my $ref (@ref){

   # parallelisation by splitting BAM file into chunks
   my $chunk = sprintf("%.0f", $ref_size{$ref} / $num_split);
   for (my $i = 0; $i < $num_split; ++$i){

      my $chunk_start = $i * $chunk;
      my $chunk_end   = ($i + 1) * $chunk;
      if ($chunk_end > $ref_size{$ref}){
         $chunk_end = $ref_size{$ref}
      }

      my $chunk = "$ref:$chunk_start-$chunk_end";
      my $str = rand_str(10);
      my $outfile = "$tmp_dir/${str}_$chunk.tsv";
      push(@chunk, $outfile);

      $manager->start and next;
      process_chunk($ref, $chunk_start, $chunk_end, $outfile);
      $manager->finish;

   }

}

$manager->wait_all_children;

warn("Merging files\n");

foreach my $infile (@chunk){
   open(IN, '<', $infile) || die "Could not open $infile: $!\n";
   while(<IN>){
      chomp;
      print "$_\n";
   }
   close(IN);
   unlink($infile);
}

warn("Done\n");

exit(0);

sub process_chunk {

   my ($ref, $chunk_start, $chunk_end, $outfile) = @_;
   my %result = ();
   my $chunk = "$ref:$chunk_start-$chunk_end";

   warn("Processing chunk: $chunk\n");

   my $sam = Bio::DB::Sam->new(
      -bam   => $bam_file,
      -fasta => $fasta_file
   );

   my @alignments = $sam->get_features_by_location(
      -seq_id => $ref,
      -start  => $chunk_start,
      -end    => $chunk_end
   );

   open(OUT, '>', $outfile) || die "Could not open $outfile for writing: $!\n";
   ALN: for my $a (@alignments) {
      my $read_id =  $a->display_name;
      if (exists $list{$read_id}){
         print OUT $a->tam_line, "\n";
      }
   }
   close(OUT);

}

sub store_read {
   my ($infile) = @_;
   my %l = ();
   open(IN, '<', $infile) || die "Could not open $infile: $!\n";
   while(<IN>){
      chomp;
      next if /^$/;
      next if /^#/;
      $l{$_} = 1;
   }
   close(IN);
   return(%l);
}

sub rand_str {
   my ($l) = @_; 
   my $s = ''; 
   my @chars = ("A".."Z", "a".."z", 0..9);
   for (1 .. $l){
      $s .= $chars[rand @chars];
   }   
   return($s)
}

sub get_ref_size {
   my ($infile) = @_;
   my %l = ();
   open(IN, '<', $infile) || die "Could not open $infile: $!\n";
   while(<IN>){
      chomp;
      next if /^\*/;
      my ($chr, $size, @rest) = split(/\t/);
      $l{$chr} = $size;
   }
   close(IN);
   return(%l);
}

sub usage {
print STDERR <<EOF;
Usage: $0 -b FILE -l FILE -t DIR -p INT -s INT

       -b infile.bam          BAM file
       -l list.txt            List of read IDs
       -f genome.fasta        FASTA file used for read alignment
       -i file.idxstats       Output from samtools idxstats saved in a file
       -t /scratch/tmp        Directory for temporary files (default /tmp/)
       -p 8                   Number of processors to use (default 8)
       -s 40                  Number of chunks to split BAM file (default 40)
       -h                     this helpful usage message

EOF
exit();
}

