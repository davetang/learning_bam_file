## README

How to extract a subset of reads (IDs stored in a file) from a BAM file?

Download BAM file with 918,571 reads.

```bash
wget -c -N http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562G1AlnRep1.bam
samtools index wgEncodeUwRepliSeqK562G1AlnRep1.bam

samtools idxstats wgEncodeUwRepliSeqK562G1AlnRep1.bam | perl -lane '$i += $F[2]; END { print $i }'
918571
```

Get some random IDs (around 0.1% of total).

```bash
samtools view -@8 -s 1984.001 wgEncodeUwRepliSeqK562G1AlnRep1.bam | cut -f1 > my_id.txt

cat my_id.txt | wc -l
     881
```

Extract sequences from `my_id.txt`.

```bash
time samtools view -@8 wgEncodeUwRepliSeqK562G1AlnRep1.bam | grep -w -f my_id.txt > my_seq.txt

real    35m48.658s
user    35m46.491s
sys     0m1.273s

cat my_seq.txt | wc -l
     881
```

Split and grep.

```bash
time for ref in $(samtools idxstats wgEncodeUwRepliSeqK562G1AlnRep1.bam | grep -v '^*' | cut -f1); do
   samtools view wgEncodeUwRepliSeqK562G1AlnRep1.bam $ref | grep -w -f my_id.txt > $ref.txt
done

real    40m12.463s
user    39m48.597s
sys     0m9.205s

cat chr*.txt | wc -l
     881
```

Split and grep using `parallel` and 7 threads.

```bash
time samtools idxstats wgEncodeUwRepliSeqK562G1AlnRep1.bam |
   grep -v '^*' |
   cut -f1 |
   parallel -j 7 --verbose 'samtools view wgEncodeUwRepliSeqK562G1AlnRep1.bam {} | grep -w -f my_id.txt > {}.txt'

real    10m26.894s
user    67m26.912s
sys     0m7.958s

cat chr*.txt | wc -l
     881
```

Back to BAM.

```bash
samtools view -H wgEncodeUwRepliSeqK562G1AlnRep1.bam > header

cat header chr*.txt | samtools view -Sb | samtools sort - > my_id.bam
```

As a Perl script using `Bio::DB::Sam` and `Parallel::ForkManager`.

```bash
samtools idxstats wgEncodeUwRepliSeqK562G1AlnRep1.bam > wgEncodeUwRepliSeqK562G1AlnRep1.idxstats

wget -c -N https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
samtools faidx hg19.fa

script/get_reads.pl -p 4 -b wgEncodeUwRepliSeqK562G1AlnRep1.bam -r 27 -l my_id.txt -i wgEncodeUwRepliSeqK562G1AlnRep1.idxstats -f hg19.fa > reads.txt
```

Use my Docker image if you have problem installing the Perl packages.

```bash
docker pull davetang/bioperl

docker run --rm -v ~/github/learning_bam_file/:/work -it davetang/bioperl /bin/bash

cd /work

perl script/get_reads.pl
Usage: script/get_reads.pl -b FILE -l FILE -t DIR -p INT -s INT -r INT

       -b infile.bam          BAM file
       -r 100                 Length of a read
       -l list.txt            List of read IDs
       -f genome.fasta        FASTA file used for read alignment
       -i file.idxstats       Output from samtools idxstats saved in a file
       -t /scratch/tmp        Directory for temporary files (default /tmp/)
       -p 8                   Number of processors to use (default 8)
       -s 40                  Number of chunks to split BAM file (default 40)
       -h                     this helpful usage message

samtools view aln.bam | cut -f1 | head -100 | sort -u > list.txt

samtools idxstats aln.bam > aln.idxstats

perl script/get_reads.pl -b aln.bam -r 100 -l list.txt -f sequence/ref.fa -i aln.idxstats -p 4 > my_reads.txt
```

