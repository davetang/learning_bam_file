#!/usr/bin/env bash

set -euo pipefail

export PATH=$PATH:miniconda3/bin/

out_md=tmp.md
miniconda3/bin/Rscript -e "rmarkdown::render('learning_bam_file.Rmd', output_file=\"$out_md\")"
github-markdown-toc/gh-md-toc $out_md > toc

cat toc <(echo) <(date) <(echo) $out_md > README.md

rm $out_md toc

>&2 echo Done!

exit 0

