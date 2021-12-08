#!/usr/bin/env bash

set -euo pipefail

out_md=tmp.md
Rscript -e "rmarkdown::render('learning_bam_file.Rmd', output_file=\"$out_md\")"
github-markdown-toc/gh-md-toc $out_md > toc

cp -f $out_md mkdocs/docs/index.md
cat toc <(echo) <(date) <(echo) $out_md > README.md

rm $out_md toc

>&2 echo Done!

exit 0

