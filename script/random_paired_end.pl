#!/usr/bin/env perl

# Simple script that takes an input fasta sequence
# and generates paired end reads

use strict;
use warnings;
use Getopt::Std;

my %opts = ();
getopts('h:f:l:n:m:s:d:', \%opts);

if ($opts{'h'} ||
    !exists $opts{'f'} ||
    !exists $opts{'l'} ||
    !exists $opts{'n'} ||
    !exists $opts{'m'}
){
   usage();
}

my $fasta = $opts{'f'};
my $len = $opts{'l'};
my $num = $opts{'n'};
my $inner_mate = $opts{'m'};
my $seed = 1984;
my $discord = 0;

if (exists $opts{'s'}){
   $seed = $opts{'s'};
}

if (exists $opts{'d'}){
   $discord = $opts{'d'};
}

srand($seed);

my $test = 0;
my $limit = 10000;
foreach my $i (1 .. $limit){
}

my $fh;
if ($fasta =~ /\.gz$/){
   open($fh, '-|', "gunzip -c $fasta") or die "Could not open $fasta $!\n";
} else {
   open($fh, '<', $fasta) or die "Could not open $fasta $!\n";
}

my %seq = ();
my $id = '';
while(<$fh>){
   chomp;
   if (/^>(.*)/){
      $id = $1;
      next;
   } else {
      $seq{$id} .= $_;
   }
}
close($fh);

my $name = 'l' . $len . '_' . 'n' . $num . '_' . 'd' . $inner_mate . '_' . $seed;
my $first_out = $name . '_1.fq.gz';
my $second_out = $name . '_2.fq.gz';

open(my $read1, '|-', "gzip >$first_out") or die "Could not write output to $first_out: $!\n";
open(my $read2, '|-', "gzip >$second_out") or die "Could not write output to  $second_out: $!\n";

for (1 .. $num){

   my @seq_id = keys %seq;
   my $seq_id = $seq_id[rand(scalar @seq_id)];
   my $seq = $seq{$seq_id};

   if (scalar @seq_id > 1){
      my $index = 0;
      $index++ until $seq_id[$index] eq $seq_id;
      splice(@seq_id, $index, 1);
   }

   my $seq_len = length($seq);
   my $limit = $seq_len - $len -$len - $inner_mate;
   if ($len > $seq_len){
      die "Your read length ($len) is longer than the sequence $seq_id\n";
   }

   # on Illumina 1.8+ ! is the worst quality
   # and J is the best
   my $fake_qual = 'J' x $len;

   my $first_start = int(rand($limit));
   my $first_read = substr($seq, $first_start, $len);
   my $first_pos = $first_start + 1;
   print $read1 "\@$_:${seq_id}_$first_pos\n$first_read\n+\n$fake_qual\n";

   my $seq_id2 = $seq_id;
   my $second_start = '';
   my $second_read = '';
   if ($discord > rand(1) && scalar @seq_id > 1){
      $seq_id2 = $seq_id[rand(scalar @seq_id)];
      my $seq2 = $seq{$seq_id2};
      my $seq_len = length($seq2);
      if ($len > $seq_len){
         die "Your read length ($len) is longer than the sequence $seq_id2\n";
      }
      my $limit = $seq_len - $len;
      $second_start = int(rand($limit));
      $second_read = substr($seq2, $second_start, $len);
   } else {
      $second_start = $first_start + $inner_mate;
      $second_read = substr($seq, $second_start, $len);
   }

   $second_read = reverse($second_read);
   $second_read =~ tr/ACGT/TGCA/;
   # reads IDs need to match!
   print $read2 "\@$_:${seq_id}_$first_pos\n$second_read\n+\n$fake_qual\n";
}

close($read1);
close($read2);

exit(0);

sub usage {
print STDERR <<EOF;
Usage: $0 -f <file.fa> -l <100> -n <100000> -m <500> [-s 1984] [-d 0]

Where:   -f         FASTA file
         -l         read lengths
         -n         number of reads
         -m         inner mate distance
         -s         seed (default: 1984)
         -d         fraction discordant (default: 0)
         -h         this helpful usage message

EOF
exit();
}

__END__

