#!/bin/bash

# make some directories
mkdir script
mkdir sequence

# clone some scripts
git clone https://gist.github.com/davetang/3e9c8ae57eded71c84de
git clone https://gist.github.com/davetang/fce12a99e514938b5725

# move scripts into script directory
mv 3e9c8ae57eded71c84de/generate_random_seq.pl script
mv fce12a99e514938b5725/random_paired_end.pl script
chmod 755 script/random_paired_end.pl

# remove gist folders
rm -rf 3e9c8ae57eded71c84de fce12a99e514938b5725

# generate random reference and reads
script/generate_random_seq.pl 1000000 31 > sequence/ref.fa
script/random_paired_end.pl sequence/ref.fa 100 10000 300 31
script/random_paired_end.pl sequence/ref.fa 100 10000 300 42
mv *.fq sequence

# get bwa
git clone https://github.com/lh3/bwa.git
cd bwa
make
cd ..

# index and align using BWA MEM
bwa/bwa index sequence/ref.fa
bwa/bwa mem sequence/ref.fa sequence/l100_n10000_d300_31_1.fq sequence/l100_n10000_d300_31_2.fq > aln.sam
bwa/bwa mem sequence/ref.fa sequence/l100_n10000_d300_42_1.fq sequence/l100_n10000_d300_42_2.fq > aln2.sam

