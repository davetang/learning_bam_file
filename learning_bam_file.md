Introduction
============

SAMtools provides various (sub)tools for manipulating alignments in the
SAM/BAM format. The SAM (Sequence Alignment/Map) format (BAM is just the
binary form of SAM) is currently the *de facto* standard for storing
large nucleotide sequence alignments. If you are dealing with
high-throughput sequencing data, at some point you will probably have to
deal with SAM/BAM files, so familiarise yourself with them! For the
latest information, please refer to the [release
notes](https://github.com/samtools/samtools/releases).

The examples below use the `ERR188273_chrX.bam` BAM file generated as
per
<a href="https://github.com/davetang/rnaseq" class="uri">https://github.com/davetang/rnaseq</a>
using the HISAT2 + StringTie2 RNA-seq pipeline.

Installing SAMtools
===================

For installing SAMtools, I recommend using the Bioconda
[samtools](https://anaconda.org/bioconda/samtools) package. I also
recommend using
[Miniconda](https://docs.conda.io/en/latest/miniconda.html) instead of
Anaconda. I wrote a [short introduction to
Conda](https://davetang.github.io/reproducible_bioinformatics/conda.html)
if you want to find learn more.

Once you have installed Miniconda, it is easy to install SAMtools.

``` bash
conda install -c bioconda samtools
```

Basic usage
===========

If you run SAMtools on the terminal without any parameters or with
`--help`, all the available utilities are listed:

``` bash
samtools --help
```

    ## 
    ## Program: samtools (Tools for alignments in the SAM format)
    ## Version: 1.6 (using htslib 1.6)
    ## 
    ## Usage:   samtools <command> [options]
    ## 
    ## Commands:
    ##   -- Indexing
    ##      dict           create a sequence dictionary file
    ##      faidx          index/extract FASTA
    ##      index          index alignment
    ## 
    ##   -- Editing
    ##      calmd          recalculate MD/NM tags and '=' bases
    ##      fixmate        fix mate information
    ##      reheader       replace BAM header
    ##      rmdup          remove PCR duplicates
    ##      targetcut      cut fosmid regions (for fosmid pool only)
    ##      addreplacerg   adds or replaces RG tags
    ##      markdup        mark duplicates
    ## 
    ##   -- File operations
    ##      collate        shuffle and group alignments by name
    ##      cat            concatenate BAMs
    ##      merge          merge sorted alignments
    ##      mpileup        multi-way pileup
    ##      sort           sort alignment file
    ##      split          splits a file by read group
    ##      quickcheck     quickly check if SAM/BAM/CRAM file appears intact
    ##      fastq          converts a BAM to a FASTQ
    ##      fasta          converts a BAM to a FASTA
    ## 
    ##   -- Statistics
    ##      bedcov         read depth per BED region
    ##      depth          compute the depth
    ##      flagstat       simple stats
    ##      idxstats       BAM index stats
    ##      phase          phase heterozygotes
    ##      stats          generate stats (former bamcheck)
    ## 
    ##   -- Viewing
    ##      flags          explain BAM flags
    ##      tview          text alignment viewer
    ##      view           SAM<->BAM<->CRAM conversion
    ##      depad          convert padded BAM to unpadded BAM

Viewing
=======

Use [bioSyntax](https://github.com/bioSyntax/bioSyntax) to prettify your
output.

``` bash
samtools view aln.bam | sam-less
```

![bioSyntax](img/sam_less.png)

Converting a SAM file to a BAM file
===================================

A BAM file is just a SAM file but stored in binary format; you should
always convert your SAM files into BAM as BAM files are smaller in size
and are faster to manipulate.

Since I don’t have a SAM file in the example folder, let’s first create
one and check out the first ten lines. Note: remember to use `-h` to
ensure the SAM file contains the sequence header information. Generally,
I recommend storing only sorted BAM files as they use less disk space
and are faster to process.

``` bash
samtools view -h eg/ERR188273_chrX.bam > eg/ERR188273_chrX.sam
```

First notice that the SAM file is much larger than the BAM file.

``` bash
ls -lh eg/ERR188273_chrX.bam eg/ERR188273_chrX.sam
```

    ## -rw-r--r--  1 dtang  staff    67M 21 Jun 20:17 eg/ERR188273_chrX.bam
    ## -rw-r--r--  1 dtang  staff   321M 28 Jun 16:41 eg/ERR188273_chrX.sam

We can use `head` to view a SAM file.

``` bash
head eg/ERR188273_chrX.sam
```

    ## @HD  VN:1.0  SO:coordinate
    ## @SQ  SN:chrX LN:156040895
    ## @PG  ID:hisat2   PN:hisat2   VN:2.2.0    CL:"/Users/dtang/github/rnaseq/hisat2/../src/hisat2-2.2.0/hisat2-align-s --wrapper basic-0 --dta -p 4 -x ../raw/chrX_data/indexes/chrX_tran -1 /tmp/4195.inpipe1 -2 /tmp/4195.inpipe2"
    ## ERR188273.4711308    73  chrX    21649   0   5S70M   =   21649   0   CGGGTGATCACGAGGTCAGGAGATCAAGACCATCCTGGCCAACACAGTGAAACCCCATCTCTACTAAAAATACAA @@@F=DDFFHGHBHIFFHIGGIFGEGHFHIGIGIFIIIGIGIGGDHIIGIIC@>DGHCHHHGHHFFFFFDEACC@ AS:i:-5 ZS:i:-5 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:70 YT:Z:UP NH:i:2
    ## ERR188273.4711308    133 chrX    21649   0   *   =   21649   0   CTACAGGTGCCCGCCACCATGCCCAGCTAATTTTTTTTGTATTTTTAGTAGAGATGGGGTTTCACTGTGTTGGCC CB@FDFFFHHGFHIJJJJIIIIIIIGGGIJGIIJJJJJJFFHIIIIGECHEHHGGHHFF?AACCDDDDDDDDBCD YT:Z:UP
    ## ERR188273.4711308    329 chrX    233717  0   5S70M   =   233717  0   CGGGTGATCACGAGGTCAGGAGATCAAGACCATCCTGGCCAACACAGTGAAACCCCATCTCTACTAAAAATACAA @@@F=DDFFHGHBHIFFHIGGIFGEGHFHIGIGIFIIIGIGIGGDHIIGIIC@>DGHCHHHGHHFFFFFDEACC@ AS:i:-5 ZS:i:-5 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:70 YT:Z:UP NH:i:2
    ## ERR188273.14904746   99  chrX    251271  60  75M =   251317  121 GAAAAATGGGCCCAGGGGACCGGCGCTCAGCATACAGAGGACCCGCGCCGGCACCTGCCTCTGAGTTCCCTTAGT @@<DDDDDFB>HHEGIIGAGIIIBGIIG@FECH<F@GIIFAE=?BCBBCBBB5@<?CBBCCCCAACDCCCCCCCC AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:75 YS:i:-2 YT:Z:CP NH:i:1
    ## ERR188273.14904746   147 chrX    251317  60  75M =   251271  -121    GCCGGCCCCTGCCTCTGAGTTCCCTTAGTACTTATTGATCATTATCGGGGAGAGGGGGATGTGGCAGGACAATAG #######B?DAHC@EGIIGGEHHGC@GFBFCEGFCIGG@EG@@H<JIEHEF@IGEHGHIIHFGHDDFDDDDD?<B AS:i:-2 ZS:i:-7 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:6A68   YS:i:0  YT:Z:CP NH:i:1
    ## ERR188273.5849805    163 chrX    265951  1   75M =   266022  146 CGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCGCCCGCCACCACGCCCGGCTAA @CCFDFFFHHHHHJJJJJJFJJJJJIJIJJJJJJJJGHIJJJJJEHIJIJGIIJJJHHFFDDEDDDDDDDDDDDD AS:i:0  ZS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:75 YS:i:0  YT:Z:CP NH:i:2
    ## ERR188273.1232356    369 chrX    265984  1   75M =   118343251   0   GAGTAGCTGGGACTACAGGCGCCCGCCACCACGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCAC @@CA@B>DC>>+@::8-755-BBBFDDEHHBGGEGHEEIJIIGIJJIGEIIIJJJIIJJIGGHHHGGFFFFF@@C AS:i:0  ZS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:75 YT:Z:UP NH:i:10

The lines starting with the “@” sign contains the header information.
The @SQ tag is the reference sequence dictionary; SN refers to the
reference sequence name and LN refers to the reference sequence length.
If you don’t see lines starting with the “@” sign, the header
information is most likely missing. If the @SQ header is absent from the
SAM file use the command below, where ref.fa is the reference fasta file
used to map the reads, to generate @SQ information.

``` bash
samtools view -bT sequence/ref.fa aln.sam > aln.bam
```

If the header information is available, we can convert a SAM file into
BAM by using `samtools view -b`. In the newer version of SAMtools the
input format is autodetected, so we no longer need the `-S` parameter.

``` bash
samtools view -b eg/ERR188273_chrX.sam > my.bam
```

Converting a BAM file to a CRAM file
====================================

Use `samtools view` with the `-T` and `-C` arguments to convert a BAM
file into CRAM.

``` bash
samtools view -T /Users/dtang/github/rnaseq/raw/chrX_data/genome/chrX.fa -C -o eg/ERR188273_chrX.cram eg/ERR188273_chrX.bam

ls -lh eg/ERR188273_chrX.[sbcr]*am
```

    ## -rw-r--r--  1 dtang  staff    67M 21 Jun 20:17 eg/ERR188273_chrX.bam
    ## -rw-r--r--  1 dtang  staff    40M 28 Jun 16:41 eg/ERR188273_chrX.cram
    ## -rw-r--r--  1 dtang  staff   321M 28 Jun 16:41 eg/ERR188273_chrX.sam

You can use `samtools view` to view a CRAM file.

``` bash
samtools view eg/ERR188273_chrX.cram | head
```

    ## ERR188273.4711308    73  chrX    21649   0   5S70M   =   21649   0   CGGGTGATCACGAGGTCAGGAGATCAAGACCATCCTGGCCAACACAGTGAAACCCCATCTCTACTAAAAATACAA @@@F=DDFFHGHBHIFFHIGGIFGEGHFHIGIGIFIIIGIGIGGDHIIGIIC@>DGHCHHHGHHFFFFFDEACC@ AS:i:-5 ZS:i:-5 XN:i:0  XM:i:0  XO:i:0  XG:i:0  YT:Z:UP NH:i:2  MD:Z:70 NM:i:0
    ## ERR188273.4711308    133 chrX    21649   0   *   =   21649   0   CTACAGGTGCCCGCCACCATGCCCAGCTAATTTTTTTTGTATTTTTAGTAGAGATGGGGTTTCACTGTGTTGGCC CB@FDFFFHHGFHIJJJJIIIIIIIGGGIJGIIJJJJJJFFHIIIIGECHEHHGGHHFF?AACCDDDDDDDDBCD YT:Z:UP
    ## ERR188273.4711308    329 chrX    233717  0   5S70M   =   233717  0   CGGGTGATCACGAGGTCAGGAGATCAAGACCATCCTGGCCAACACAGTGAAACCCCATCTCTACTAAAAATACAA @@@F=DDFFHGHBHIFFHIGGIFGEGHFHIGIGIFIIIGIGIGGDHIIGIIC@>DGHCHHHGHHFFFFFDEACC@ AS:i:-5 ZS:i:-5 XN:i:0  XM:i:0  XO:i:0  XG:i:0  YT:Z:UP NH:i:2  MD:Z:70 NM:i:0
    ## ERR188273.14904746   99  chrX    251271  60  75M =   251317  121 GAAAAATGGGCCCAGGGGACCGGCGCTCAGCATACAGAGGACCCGCGCCGGCACCTGCCTCTGAGTTCCCTTAGT @@<DDDDDFB>HHEGIIGAGIIIBGIIG@FECH<F@GIIFAE=?BCBBCBBB5@<?CBBCCCCAACDCCCCCCCC AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  YS:i:-2 YT:Z:CP NH:i:1  MD:Z:75 NM:i:0
    ## ERR188273.14904746   147 chrX    251317  60  75M =   251271  -121    GCCGGCCCCTGCCTCTGAGTTCCCTTAGTACTTATTGATCATTATCGGGGAGAGGGGGATGTGGCAGGACAATAG #######B?DAHC@EGIIGGEHHGC@GFBFCEGFCIGG@EG@@H<JIEHEF@IGEHGHIIHFGHDDFDDDDD?<B AS:i:-2 ZS:i:-7 XN:i:0  XM:i:1  XO:i:0  XG:i:0  YS:i:0  YT:Z:CP NH:i:1  MD:Z:6A68   NM:i:1
    ## ERR188273.5849805    163 chrX    265951  1   75M =   266022  146 CGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCGCCCGCCACCACGCCCGGCTAA @CCFDFFFHHHHHJJJJJJFJJJJJIJIJJJJJJJJGHIJJJJJEHIJIJGIIJJJHHFFDDEDDDDDDDDDDDD AS:i:0  ZS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  YS:i:0  YT:Z:CP NH:i:2  MD:Z:75 NM:i:0
    ## ERR188273.1232356    369 chrX    265984  1   75M =   118343251   0   GAGTAGCTGGGACTACAGGCGCCCGCCACCACGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCAC @@CA@B>DC>>+@::8-755-BBBFDDEHHBGGEGHEEIJIIGIJJIGEIIIJJJIIJJIGGHHHGGFFFFF@@C AS:i:0  ZS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  YT:Z:UP NH:i:10 MD:Z:75 NM:i:0
    ## ERR188273.5927795    385 chrX    265991  1   75M =   114048277   0   TGGGACTACAGGCGCCCGCCACCACGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTA =?BB??BD?FBHHBEAE@CDGG@HH=FA@GEGE;FGACCHBE6?A=ACE9)7@DCE>>5'3=338:;:>2<AA?: AS:i:0  ZS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  YT:Z:UP NH:i:10 MD:Z:75 NM:i:0
    ## ERR188273.5849805    83  chrX    266022  1   75M =   265951  -146    CTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTGTCGATCTCCTGACCTCGTG DDDDDDDEEEEEDBFEDGHHHHHHJHIJJIGIGHFBJJIHGJJIIJJJJJJJJIGIJJJJJJHHHHHDFFFFCCB AS:i:0  ZS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  YS:i:0  YT:Z:CP NH:i:2  MD:Z:75 NM:i:0
    ## ERR188273.13655123   113 chrX    266022  1   75M =   118343234   0   CTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTGTCGATCTCCTGACCTCGTG AACBBBACCCC>;3?BCFFEEHHHEEGIGGHAGFBBHFBHHEHCG@<@ABG??@@?BB9GBGAFFD<<DDAD@@@ AS:i:0  ZS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  YT:Z:UP NH:i:2  MD:Z:75 NM:i:0

I have an [old blog
post](https://davetang.org/muse/2014/09/26/bam-to-cram/) on the CRAM
format.

Sorting a SAM/BAM file
======================

Always sort your SAM/BAM files; many downstream programs only take
sorted BAM files. In SAMtools version 1.3 or newer, you can sort a SAM
file directly.

``` bash
samtools sort eg/ERR188273_chrX.sam -o eg/sorted.bam
```

You should use use additional threads by specifying `-@ 4` (using 4
threads) to speed up sorting.

``` bash
time samtools sort eg/ERR188273_chrX.sam -o eg/sorted.bam
time samtools sort -@ 4 eg/ERR188273_chrX.sam -o eg/sorted.bam
```

    ## 
    ## real 0m12.119s
    ## user 0m11.654s
    ## sys  0m0.368s
    ## [bam_sort_core] merging from 0 files and 4 in-memory blocks...
    ## 
    ## real 0m4.160s
    ## user 0m12.334s
    ## sys  0m0.442s

Creating a BAM index file
=========================

Various tools require BAM index files, such as IGV, which is a program
for visualising a BAM file.

``` bash
samtools index eg/ERR188273_chrX.bam
```

Filtering unmapped reads
========================

Use `-F 4` to filter out unmapped reads. Use the `flags` subcommand to
find out what a flag represents.

``` bash
samtools flags 4
samtools view -F 4 -b eg/ERR188273_chrX.bam > eg/ERR188273_chrX.mapped.bam
```

    ## 0x4  4   UNMAP

Use `-f 4` to keep only unmapped reads.

``` bash
samtools view -f 4 -b eg/ERR188273_chrX.bam > eg/ERR188273_chrX.unmapped.bam
```

Extracting entries mapping to a specific loci
=============================================

If we want all reads mapping within a specific genomic region, we can
use `samtools view` and the `ref:start-end` syntax. You can use just the
`ref` to extract an entire reference sequence such as a chromosome
(example not shown here). This requires a BAM index file.

``` bash
samtools view eg/ERR188273_chrX.bam chrX:20000-30000
```

    ## ERR188273.4711308    73  chrX    21649   0   5S70M   =   21649   0   CGGGTGATCACGAGGTCAGGAGATCAAGACCATCCTGGCCAACACAGTGAAACCCCATCTCTACTAAAAATACAA @@@F=DDFFHGHBHIFFHIGGIFGEGHFHIGIGIFIIIGIGIGGDHIIGIIC@>DGHCHHHGHHFFFFFDEACC@ AS:i:-5 ZS:i:-5 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:70 YT:Z:UP NH:i:2
    ## ERR188273.4711308    133 chrX    21649   0   *   =   21649   0   CTACAGGTGCCCGCCACCATGCCCAGCTAATTTTTTTTGTATTTTTAGTAGAGATGGGGTTTCACTGTGTTGGCC CB@FDFFFHHGFHIJJJJIIIIIIIGGGIJGIIJJJJJJFFHIIIIGECHEHHGGHHFF?AACCDDDDDDDDBCD YT:Z:UP

Note that this takes into account the mapping of the entire read and not
just the starting position. For example, if you specified
chrX:20000-30000, a read that is 75 bp long that maps to position 19999
will also be returned. You can save the output as another BAM file if
you wish.

``` bash
samtools view -b eg/ERR188273_chrX.bam chrX:20000-30000 > eg/ERR188273_chrX_20000_30000.bam
```

You can also use a BED file, with several entries, to extract reads of
interest.

``` bash
cat eg/my.bed 

samtools view -L eg/my.bed eg/ERR188273_chrX.bam
```

    ## chrX 20000   30000
    ## chrX 233000  260000
    ## ERR188273.4711308    73  chrX    21649   0   5S70M   =   21649   0   CGGGTGATCACGAGGTCAGGAGATCAAGACCATCCTGGCCAACACAGTGAAACCCCATCTCTACTAAAAATACAA @@@F=DDFFHGHBHIFFHIGGIFGEGHFHIGIGIFIIIGIGIGGDHIIGIIC@>DGHCHHHGHHFFFFFDEACC@ AS:i:-5 ZS:i:-5 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:70 YT:Z:UP NH:i:2
    ## ERR188273.4711308    133 chrX    21649   0   *   =   21649   0   CTACAGGTGCCCGCCACCATGCCCAGCTAATTTTTTTTGTATTTTTAGTAGAGATGGGGTTTCACTGTGTTGGCC CB@FDFFFHHGFHIJJJJIIIIIIIGGGIJGIIJJJJJJFFHIIIIGECHEHHGGHHFF?AACCDDDDDDDDBCD YT:Z:UP
    ## ERR188273.4711308    329 chrX    233717  0   5S70M   =   233717  0   CGGGTGATCACGAGGTCAGGAGATCAAGACCATCCTGGCCAACACAGTGAAACCCCATCTCTACTAAAAATACAA @@@F=DDFFHGHBHIFFHIGGIFGEGHFHIGIGIFIIIGIGIGGDHIIGIIC@>DGHCHHHGHHFFFFFDEACC@ AS:i:-5 ZS:i:-5 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:70 YT:Z:UP NH:i:2
    ## ERR188273.14904746   99  chrX    251271  60  75M =   251317  121 GAAAAATGGGCCCAGGGGACCGGCGCTCAGCATACAGAGGACCCGCGCCGGCACCTGCCTCTGAGTTCCCTTAGT @@<DDDDDFB>HHEGIIGAGIIIBGIIG@FECH<F@GIIFAE=?BCBBCBBB5@<?CBBCCCCAACDCCCCCCCC AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:75 YS:i:-2 YT:Z:CP NH:i:1
    ## ERR188273.14904746   147 chrX    251317  60  75M =   251271  -121    GCCGGCCCCTGCCTCTGAGTTCCCTTAGTACTTATTGATCATTATCGGGGAGAGGGGGATGTGGCAGGACAATAG #######B?DAHC@EGIIGGEHHGC@GFBFCEGFCIGG@EG@@H<JIEHEF@IGEHGHIIHFGHDDFDDDDD?<B AS:i:-2 ZS:i:-7 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:6A68   YS:i:0  YT:Z:CP NH:i:1

Extracting only the first read from paired end BAM files
========================================================

Sometimes you only want the first pair of a mate. 0x0040 is hexadecimal
for 64 (i.e. 16 \* 4), which is binary for 1000000, corresponding to the
read in the first read pair.

``` bash
samtools view -b -f 0x0040 eg/ERR188273_chrX.bam > eg/first.bam
```

Stats
=====

For simple statistics use `samtools flagstat`.

``` bash
samtools flagstat eg/ERR188273_chrX.bam
```

    ## 1176360 + 0 in total (QC-passed reads + QC-failed reads)
    ## 16276 + 0 secondary
    ## 0 + 0 supplementary
    ## 0 + 0 duplicates
    ## 1126961 + 0 mapped (95.80% : N/A)
    ## 1160084 + 0 paired in sequencing
    ## 580042 + 0 read1
    ## 580042 + 0 read2
    ## 1060858 + 0 properly paired (91.45% : N/A)
    ## 1065618 + 0 with itself and mate mapped
    ## 45067 + 0 singletons (3.88% : N/A)
    ## 0 + 0 with mate mapped to a different chr
    ## 0 + 0 with mate mapped to a different chr (mapQ>=5)

For additional stats, use `samtools stats`.

``` bash
samtools stats eg/ERR188273_chrX.bam | grep ^SN
```

    ## SN   raw total sequences:    1160084
    ## SN   filtered sequences: 0
    ## SN   sequences:  1160084
    ## SN   is sorted:  1
    ## SN   1st fragments:  580042
    ## SN   last fragments: 580042
    ## SN   reads mapped:   1110685
    ## SN   reads mapped and paired:    1065618 # paired-end technology bit set + both mates mapped
    ## SN   reads unmapped: 49399
    ## SN   reads properly paired:  1060858 # proper-pair bit set
    ## SN   reads paired:   1160084 # paired-end technology bit set
    ## SN   reads duplicated:   0   # PCR or optical duplicate bit set
    ## SN   reads MQ0:  905 # mapped and MQ=0
    ## SN   reads QC failed:    0
    ## SN   non-primary alignments: 16276
    ## SN   total length:   87006300    # ignores clipping
    ## SN   bases mapped:   83301375    # ignores clipping
    ## SN   bases mapped (cigar):   83064942    # more accurate
    ## SN   bases trimmed:  0
    ## SN   bases duplicated:   0
    ## SN   mismatches: 423271  # from NM fields
    ## SN   error rate: 5.095663e-03    # mismatches / bases mapped (cigar)
    ## SN   average length: 75
    ## SN   maximum length: 75
    ## SN   average quality:    36.0
    ## SN   insert size average:    182.6
    ## SN   insert size standard deviation: 175.5
    ## SN   inward oriented pairs:  530549
    ## SN   outward oriented pairs: 834
    ## SN   pairs with other orientation:   1004
    ## SN   pairs on different chromosomes: 0

Interpreting the BAM flags
==========================

The second column in a SAM/BAM file is the flag column. They may seem
confusing at first but the encoding allows details about a read to be
stored by just using a few digits. The trick is to convert the numerical
digit into binary, and then use the table to interpret the binary
numbers, where 1 = true and 0 = false. I wrote a blog post on BAM flags:
<a href="http://davetang.org/muse/2014/03/06/understanding-bam-flags/" class="uri">http://davetang.org/muse/2014/03/06/understanding-bam-flags/</a>,
which also includes a Perl script for interpreting BAM flags. There is
also the `flags` subcommand.

``` bash
samtools flags
```

    ## 
    ## About: Convert between textual and numeric flag representation
    ## Usage: samtools flags INT|STR[,...]
    ## 
    ## Flags:
    ##  0x1 PAIRED        .. paired-end (or multiple-segment) sequencing technology
    ##  0x2 PROPER_PAIR   .. each segment properly aligned according to the aligner
    ##  0x4 UNMAP         .. segment unmapped
    ##  0x8 MUNMAP        .. next segment in the template unmapped
    ##  0x10    REVERSE       .. SEQ is reverse complemented
    ##  0x20    MREVERSE      .. SEQ of the next segment in the template is reversed
    ##  0x40    READ1         .. the first segment in the template
    ##  0x80    READ2         .. the last segment in the template
    ##  0x100   SECONDARY     .. secondary alignment
    ##  0x200   QCFAIL        .. not passing quality controls
    ##  0x400   DUP           .. PCR or optical duplicate
    ##  0x800   SUPPLEMENTARY .. supplementary alignment

samtools calmd/fillmd
=====================

The calmd or fillmd tool is useful for visualising mismatches and
insertions in an alignment of a read to a reference genome. The `-e`
argument changes identical bases between the read and reference into
`=`.

``` bash
samtools view -b eg/ERR188273_chrX.bam | samtools fillmd -e - /Users/dtang/github/rnaseq/raw/chrX_data/genome/chrX.fa > eg/ERR188273_chrX_fillmd.bam

head eg/ERR188273_chrX_fillmd.bam
```

    ## @HD  VN:1.0  SO:coordinate
    ## @SQ  SN:chrX LN:156040895
    ## @PG  ID:hisat2   PN:hisat2   VN:2.2.0    CL:"/Users/dtang/github/rnaseq/hisat2/../src/hisat2-2.2.0/hisat2-align-s --wrapper basic-0 --dta -p 4 -x ../raw/chrX_data/indexes/chrX_tran -1 /tmp/4195.inpipe1 -2 /tmp/4195.inpipe2"
    ## ERR188273.4711308    73  chrX    21649   0   5S70M   =   21649   0   CGGGT====================================================================== @@@F=DDFFHGHBHIFFHIGGIFGEGHFHIGIGIFIIIGIGIGGDHIIGIIC@>DGHCHHHGHHFFFFFDEACC@ AS:i:-5 ZS:i:-5 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:70 YT:Z:UP NH:i:2
    ## ERR188273.4711308    133 chrX    21649   0   *   =   21649   0   CTACAGGTGCCCGCCACCATGCCCAGCTAATTTTTTTTGTATTTTTAGTAGAGATGGGGTTTCACTGTGTTGGCC CB@FDFFFHHGFHIJJJJIIIIIIIGGGIJGIIJJJJJJFFHIIIIGECHEHHGGHHFF?AACCDDDDDDDDBCD YT:Z:UP
    ## ERR188273.4711308    329 chrX    233717  0   5S70M   =   233717  0   CGGGT====================================================================== @@@F=DDFFHGHBHIFFHIGGIFGEGHFHIGIGIFIIIGIGIGGDHIIGIIC@>DGHCHHHGHHFFFFFDEACC@ AS:i:-5 ZS:i:-5 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:70 YT:Z:UP NH:i:2
    ## ERR188273.14904746   99  chrX    251271  60  75M =   251317  121 =========================================================================== @@<DDDDDFB>HHEGIIGAGIIIBGIIG@FECH<F@GIIFAE=?BCBBCBBB5@<?CBBCCCCAACDCCCCCCCC AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:75 YS:i:-2 YT:Z:CP NH:i:1
    ## ERR188273.14904746   147 chrX    251317  60  75M =   251271  -121    ======C==================================================================== #######B?DAHC@EGIIGGEHHGC@GFBFCEGFCIGG@EG@@H<JIEHEF@IGEHGHIIHFGHDDFDDDDD?<B AS:i:-2 ZS:i:-7 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:6A68   YS:i:0  YT:Z:CP NH:i:1
    ## ERR188273.5849805    163 chrX    265951  1   75M =   266022  146 =========================================================================== @CCFDFFFHHHHHJJJJJJFJJJJJIJIJJJJJJJJGHIJJJJJEHIJIJGIIJJJHHFFDDEDDDDDDDDDDDD AS:i:0  ZS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:75 YS:i:0  YT:Z:CP NH:i:2
    ## ERR188273.1232356    369 chrX    265984  1   75M =   118343251   0   =========================================================================== @@CA@B>DC>>+@::8-755-BBBFDDEHHBGGEGHEEIJIIGIJJIGEIIIJJJIIJJIGGHHHGGFFFFF@@C AS:i:0  ZS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:75 YT:Z:UP NH:i:10

Creating fastq files from a BAM file
====================================

``` bash
samtools fastq -1 eg/ERR188273_chrX_1.fq -2 eg/ERR188273_chrX_2.fq eg/ERR188273_chrX.bam
head eg/ERR188273_chrX_1.fq
```

    ## [M::bam2fq_mainloop] discarded 0 singletons
    ## [M::bam2fq_mainloop] processed 1160084 reads
    ## @ERR188273.4711308
    ## CGGGTGATCACGAGGTCAGGAGATCAAGACCATCCTGGCCAACACAGTGAAACCCCATCTCTACTAAAAATACAA
    ## +
    ## @@@F=DDFFHGHBHIFFHIGGIFGEGHFHIGIGIFIIIGIGIGGDHIIGIIC@>DGHCHHHGHHFFFFFDEACC@
    ## @ERR188273.14904746
    ## GAAAAATGGGCCCAGGGGACCGGCGCTCAGCATACAGAGGACCCGCGCCGGCACCTGCCTCTGAGTTCCCTTAGT
    ## +
    ## @@<DDDDDFB>HHEGIIGAGIIIBGIIG@FECH<F@GIIFAE=?BCBBCBBB5@<?CBBCCCCAACDCCCCCCCC
    ## @ERR188273.5849805
    ## CACGAGGTCAGGAGATCGACACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAG

Random subsampling of BAM file
==============================

The SAMtools view `-s` parameter allows you to randomly sample lines of
a BAM file. Using `0.5` will subsample half of all mapped reads.

``` bash
samtools view -s 0.5 -b eg/ERR188273_chrX.bam > eg/ERR188273_chrX_rand.bam
```

Count number of reads
=====================

Use `samtools idxstats` to print stats on a BAM file; this requires an
index file which is created by running `samtools index`.

``` bash
# output of idxstats is:
# ref name, sequence length of ref, no. mapped reads, and no. unmapped reads
samtools idxstats eg/ERR188273_chrX.bam
```

    ## chrX 156040895   1126961 45067
    ## *    0   0   4332

We can use this with `awk` to sum up the columns.

``` bash
# number of reads = mapped + unmapped
samtools idxstats eg/ERR188273_chrX.bam | awk '{s+=$3+$4} END {print s}'

# number of mapped reads = 3rd column
samtools idxstats eg/ERR188273_chrX.bam  | awk '{s+=$3} END {print s}'
```

    ## 1176360
    ## 1126961

Obtaining genomic sequence
==========================

Use `faidx` to fetch genomic sequence; coordinates are 1-based.

``` bash
# index fasta file
samtools faidx /Users/dtang/github/rnaseq/raw/chrX_data/genome/chrX.fa

# obtain sequence
samtools faidx /Users/dtang/github/rnaseq/raw/chrX_data/genome/chrX.fa chrX:300000-300100
```

    ## >chrX:300000-300100
    ## ctgagatcgtgccactgcactccagcctgggcgacagagcgagactccatctcaaaaaaa
    ## aaaaaaaaaaaaaagaTggggtctctctatgttggccaggt

Comparing BAM files
===================

Install [deepTools](https://deeptools.readthedocs.io/en/develop/) and
use
[bamCompare](https://deeptools.readthedocs.io/en/develop/content/tools/bamCompare.html).
The bigWig output file shows the ratio of reads between `b1` and `b2` in
50 bp (default) windows.

Converting reference names
==========================

One of the most annoying bioinformatics problems is the use of different
chromosome names, e.g. chr1 vs 1, in different references even when the
sequences are identical. The GRCh38 reference downloaded from Ensembl
has chromosome names without the `chr`:

    >1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF

Whereas the reference names from UCSC has the `chr`:

    >chr1  AC:CM000663.2  gi:568336023  LN:248956422  rl:Chromosome  M5:6aef897c3d6ff0c78aff06ac189178dd  AS:GRCh38

Luckily you can change the reference names using `samtools reheader` but
just make sure your reference sequences are actually identical.

``` bash
samtools view eg/ERR188273_chrX.bam | head -2

# view header
samtools view -H eg/ERR188273_chrX.bam

# substitute header with new name
samtools view -H eg/ERR188273_chrX.bam | sed 's/SN:chrX/SN:X/' > eg/my_header

# save bam file with new ref
samtools reheader eg/my_header eg/ERR188273_chrX.bam > eg/ERR188273_X.bam

samtools view eg/ERR188273_X.bam | head -2
```

    ## ERR188273.4711308    73  chrX    21649   0   5S70M   =   21649   0   CGGGTGATCACGAGGTCAGGAGATCAAGACCATCCTGGCCAACACAGTGAAACCCCATCTCTACTAAAAATACAA @@@F=DDFFHGHBHIFFHIGGIFGEGHFHIGIGIFIIIGIGIGGDHIIGIIC@>DGHCHHHGHHFFFFFDEACC@ AS:i:-5 ZS:i:-5 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:70 YT:Z:UP NH:i:2
    ## ERR188273.4711308    133 chrX    21649   0   *   =   21649   0   CTACAGGTGCCCGCCACCATGCCCAGCTAATTTTTTTTGTATTTTTAGTAGAGATGGGGTTTCACTGTGTTGGCC CB@FDFFFHHGFHIJJJJIIIIIIIGGGIJGIIJJJJJJFFHIIIIGECHEHHGGHHFF?AACCDDDDDDDDBCD YT:Z:UP
    ## @HD  VN:1.0  SO:coordinate
    ## @SQ  SN:chrX LN:156040895
    ## @PG  ID:hisat2   PN:hisat2   VN:2.2.0    CL:"/Users/dtang/github/rnaseq/hisat2/../src/hisat2-2.2.0/hisat2-align-s --wrapper basic-0 --dta -p 4 -x ../raw/chrX_data/indexes/chrX_tran -1 /tmp/4195.inpipe1 -2 /tmp/4195.inpipe2"
    ## ERR188273.4711308    73  X   21649   0   5S70M   =   21649   0   CGGGTGATCACGAGGTCAGGAGATCAAGACCATCCTGGCCAACACAGTGAAACCCCATCTCTACTAAAAATACAA @@@F=DDFFHGHBHIFFHIGGIFGEGHFHIGIGIFIIIGIGIGGDHIIGIIC@>DGHCHHHGHHFFFFFDEACC@ AS:i:-5 ZS:i:-5 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:70 YT:Z:UP NH:i:2
    ## ERR188273.4711308    133 X   21649   0   *   =   21649   0   CTACAGGTGCCCGCCACCATGCCCAGCTAATTTTTTTTGTATTTTTAGTAGAGATGGGGTTTCACTGTGTTGGCC CB@FDFFFHHGFHIJJJJIIIIIIIGGGIJGIIJJJJJJFFHIIIIGECHEHHGGHHFF?AACCDDDDDDDDBCD YT:Z:UP

Coverage
========

We can use `samtools depth` to tally the number of reads covering a
region; the three columns are the reference, position, and read
coverage. In the example below, there are two reads covering positions
200 - 205. The `samtools mpileup` command can provide more information,
including:

1.  Sequence name
2.  1-based coordinate
3.  Reference base
4.  Number of reads covering this position
5.  Read bases
6.  Base qualities
7.  Alignment mapping qualities

See
<a href="https://davetang.org/muse/2015/08/26/samtools-mpileup/" class="uri">https://davetang.org/muse/2015/08/26/samtools-mpileup/</a>
for more information.

``` bash
samtools depth eg/ERR188273_chrX.bam | head
```

    ## chrX 21649   1
    ## chrX 21650   1
    ## chrX 21651   1
    ## chrX 21652   1
    ## chrX 21653   1
    ## chrX 21654   1
    ## chrX 21655   1
    ## chrX 21656   1
    ## chrX 21657   1
    ## chrX 21658   1

Stargazers over time
====================

[![Stargazers over
time](https://starchart.cc/davetang/learning_bam_file.svg)](https://starchart.cc/davetang/learning_bam_file)
