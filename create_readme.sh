#!/usr/bin/env bash

if ! [ -x "$(command -v gh-md-toc)" ]; then
  >&2 echo 'Please install gh-md-toc from https://github.com/ekalinin/github-markdown-toc'
  exit 1
fi

if ! [ -x "$(command -v Rscript)" ]; then
  >&2 echo 'Please make sure you have R installed'
  exit 1
fi

if ! [ -f ~/github/rnaseq/raw/chrX_data/genome/chrX.fa ]; then
  >&2 echo 'Please make sure you have cloned https://github.com/davetang/rnaseq and downloaded the reference files'
  exit 1
fi

out_md=tmp.md
Rscript -e "rmarkdown::render('learning_bam_file.Rmd', output_file=\"$out_md\")"
gh-md-toc $out_md > toc

cat toc <(echo) <(date) <(echo) $out_md > README.md

rm $out_md toc

echo Done!

exit 0

