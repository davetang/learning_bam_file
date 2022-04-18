#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "Usage: $0 <num> <bp> <seed>\n";
my $num = shift or die $usage;
my $len = shift or die $usage;
my $seed = shift or die $usage;

# set seed for reproducibility
srand($seed);

foreach my $i (1 .. $num){
   my $random_seq = random_seq($len);
   print ">${i}\n$random_seq\n";
}

exit(0);

sub random_seq {
   my ($len) = @_;
   my @nuc = qw/ A C G T /;
   my $seq = '';
   for (1 .. $len){
      my $rand_ind = int(rand(scalar(@nuc)));
      $seq .= $nuc[$rand_ind];
   }
   return($seq);
}

exit(0);
