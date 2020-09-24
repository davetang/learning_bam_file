## README

[bam-readcount](https://github.com/genome/bam-readcount) tool for generating metrics at single nucleotide positions from BAM files. We will use Docker to install `bam-readcount`; the Docker image is available at https://hub.docker.com/r/davetang/base. If you have never used Docker before, I wrote a [short guide](https://davetang.github.io/reproducible_bioinformatics/docker.html).

```bash
cd ${HOME}/github/learning_bam_file

# image available at https://hub.docker.com/r/davetang/base
docker run --rm -it -v $(pwd):/data davetang/base /bin/bash

# run the following commands inside the Docker container
mkdir ~/github/ && cd ~/github/
git clone --recursive https://github.com/genome/bam-readcount.git

# build and compile
mkdir -p ~/tool/ && cd ~/tool/
cmake ~/github/bam-readcount
make

# check version
bin/bam-readcount --version
bam-readcount version: 0.8.0-unstable-7-625eea2 (commit 625eea2)

# check usage
bin/bam-readcount 
Usage: bam-readcount [OPTIONS] <bam_file> [region]
Generate metrics for bam_file at single nucleotide positions.
Example: bam-readcount -f ref.fa some.bam

Available options:
  -h [ --help ]                         produce this message
  -v [ --version ]                      output the version number
  -q [ --min-mapping-quality ] arg (=0) minimum mapping quality of reads used 
                                        for counting.
  -b [ --min-base-quality ] arg (=0)    minimum base quality at a position to 
                                        use the read for counting.
  -d [ --max-count ] arg (=10000000)    max depth to avoid excessive memory 
                                        usage.
  -l [ --site-list ] arg                file containing a list of regions to 
                                        report readcounts within.
  -f [ --reference-fasta ] arg          reference sequence in the fasta format.
  -D [ --print-individual-mapq ] arg    report the mapping qualities as a comma
                                        separated list.
  -p [ --per-library ]                  report results by library.
  -w [ --max-warnings ] arg             maximum number of warnings of each type
                                        to emit. -1 gives an unlimited number.
  -i [ --insertion-centric ]            generate indel centric readcounts. 
                                        Reads containing insertions will not be
                                        included in per-base counts
```

Prepare reference file and run `bam-readcount` on example files provided in the repo.

```bash
cd /data/eg
cp chrX.fa.bz2 ref.fa.bz2
bunzip2 ref.fa.bz2

~/tool/bin/bam-readcount -f ref.fa -w 0 -l ../bam-readcount/region.bed ERR188273_chrX.bam > ../bam-readcount/ERR188273_chrX.metrics.txt

# clean up
rm ref.fa*
```

[Output format](https://github.com/genome/bam-readcount#normal-output):

1. chr
2. position
3. reference_base
4. depth
5. base:count:avg_mapping_quality:avg_basequality:avg_se_mapping_quality:num_plus_strand:num_minus_strand:avg_pos_as_fraction:avg_num_mismatches_as_fraction:avg_sum_mismatch_qualities:num_q2_containing_reads:avg_distance_to_q2_start_in_q2_reads:avg_clipped_length:avg_distance_to_effective_3p_end

```bash
chrX    21649   g       1       =:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  A:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  C:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  G:1:0.00:35.00:0.00:1:0:0.00:0.00:0.00:1:0.91:70.00:0.91      T:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00
chrX    21650   a       1       =:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  A:1:0.00:35.00:0.00:1:0:0.03:0.00:0.00:1:0.89:70.00:0.89        C:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  G:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00        T:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00
chrX    21651   t       1       =:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  A:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  C:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  G:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00        T:1:0.00:37.00:0.00:1:0:0.06:0.00:0.00:1:0.88:70.00:0.88        N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00
```

