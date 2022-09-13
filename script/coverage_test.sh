#!/usr/bin/env bash

set -euo pipefail

bin=$(dirname $0)
root=${bin}/..

if [[ ! -d ${root}/test ]]; then
   mkdir ${root}/test
fi

outdir=${root}/test

check_depend (){
   tool=$1
   if [[ ! -x $(command -v ${tool}) ]]; then
      >&2 echo Could not find ${tool}
      exit 1
   fi
}

for t in samtools bwa column; do
   check_depend ${t}
done

ref_num=2
ref_size=1000000
seed=1984
ref=ref_${ref_num}_${ref_size}_${seed}.fa
base=$(basename ${ref} .fa)

if [[ ! -e ${outdir}/${ref} ]]; then
   >&2 echo Generating random reference
   ${bin}/generate_random_seq.pl ${ref_num} ${ref_size} ${seed} > ${outdir}/${ref}
   >&2 echo Indexing ${ref}
   bwa index ${outdir}/${ref}
fi

len=150
num=111111
md=500
read1=l${len}_n${num}_d${md}_${seed}_1.fq.gz
read2=l${len}_n${num}_d${md}_${seed}_2.fq.gz

if [[ ! -e ${outdir}/${read1} ]]; then
   >&2 echo Generating random reads
   ${bin}/random_paired_end.pl \
      -f ${outdir}/ref.fa \
      -l ${len} \
      -n ${num} \
      -m ${md} \
      -s ${seed}
   mv -f ${read1} ${read2} ${outdir}
fi

if [[ ! -e ${outdir}/${base}.bwa.bam && ! -e ${outdir}/${base}.bwa.bam.bai ]]; then
   >&2 echo Mapping reads
   bwa mem ${outdir}/${ref} ${outdir}/${read1} ${outdir}/${read2} |\
      samtools sort -O BAM |\
      tee ${outdir}/${base}.bwa.bam |\
      samtools index - ${outdir}/${base}.bwa.bam.bai
fi

cov=$(bc -l<<<"${len}*${num}*${ref_num}/(${ref_size}*${ref_num})")
>&2 echo -e "Coverage should be ${cov}\n"
>&2 echo Coverage calculation using samtools coverage
samtools coverage ${outdir}/${base}.bwa.bam | column -t
>&2 echo -e

>&2 echo Coverage calculation using samtools depth
samtools depth ${outdir}/${base}.bwa.bam | perl -ane '$t += $F[2]; END {$cov = $t / $.; printf "Bases covered:\t%.2f\nCoverage:\t%.2f\n", $., $cov}'
>&2 echo -e

>&2 echo Coverage calculation using samtools mpileup
samtools mpileup ${outdir}/${base}.bwa.bam | perl -ane '$t += $F[3]; END {$cov = $t / $.; printf "Bases covered:\t%.2f\nCoverage:\t%.2f\n", $., $cov}'
>&2 echo -e

>&2 echo Done
exit 0

