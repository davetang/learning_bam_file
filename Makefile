all: readme
data: genome/chrX.fa
tools: github-markdown-toc samtools pandoc
pandoc: /usr/bin/pandoc
samtools_ver = 1.14

genome/chrX.fa:
	bunzip2 -c eg/chrX.fa.bz2 > genome/chrX.fa

github-markdown-toc:
	git clone https://github.com/ekalinin/github-markdown-toc.git

samtools:
	wget https://github.com/samtools/samtools/releases/download/$(samtools_ver)/samtools-$(samtools_ver).tar.bz2 && tar xjf samtools-$(samtools_ver).tar.bz2 && cd samtools-$(samtools_ver) && ./configure && make && mv samtools .. && cd .. && rm -rf samtools-$(samtools_ver).tar.bz2 samtools-$(samtools_ver)

/usr/bin/pandoc:
	apt update && apt install -y pandoc

readme: data tools
	./create_readme.sh

clean:
	rm -rf genome/chrX.fa github-markdown-toc samtools tmp.html && eg/clean.sh
