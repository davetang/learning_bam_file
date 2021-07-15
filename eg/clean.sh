#!/usr/bin/env bash

set -euo pipefail

dir=$(dirname $0)

rm -f \
  ${dir}/ERR188273_X.bam \
  ${dir}/ERR188273_chrX.cram \
  ${dir}/ERR188273_chrX.mapped.bam \
  ${dir}/ERR188273_chrX.sam \
  ${dir}/ERR188273_chrX.stats \
  ${dir}/ERR188273_chrX.unmapped.bam \
  ${dir}/ERR188273_chrX_1.fq \
  ${dir}/ERR188273_chrX_2.fq \
  ${dir}/ERR188273_chrX_20000_30000.bam \
  ${dir}/ERR188273_chrX_fillmd.bam \
  ${dir}/ERR188273_chrX_rand.bam \
  ${dir}/first.bam \
  ${dir}/my.bam \
  ${dir}/my_header \
  ${dir}/sorted.bam

