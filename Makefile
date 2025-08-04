.PHONY: all data tools pandoc clean

all: readme
data: genome/chrX.fa
tools: github-markdown-toc samtools bwa minimap2 mosdepth pandoc
pandoc: /usr/bin/pandoc
samtools_ver := 1.22.1
minimap2_ver := 2.24
mosdepth_ver := 0.3.7

genome/chrX.fa:
	bunzip2 -c eg/chrX.fa.bz2 > $@

github-markdown-toc:
	git clone https://github.com/ekalinin/github-markdown-toc.git

samtools:
	wget --quiet https://github.com/samtools/samtools/releases/download/$(samtools_ver)/samtools-$(samtools_ver).tar.bz2 \
		&& tar xjf samtools-$(samtools_ver).tar.bz2 \
		&& cd samtools-$(samtools_ver) \
		&& ./configure \
		&& make \
		&& mv samtools .. \
		&& cd .. \
		&& rm -rf samtools-$(samtools_ver).tar.bz2 samtools-$(samtools_ver)

bwa:
	wget --quiet https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 \
		&& tar xjf bwa-0.7.17.tar.bz2 \
		&& cd bwa-0.7.17 \
		&& make \
		&& mv bwa .. \
		&& cd .. \
		&& rm -rf bwa-*

minimap2:
	wget --quiet https://github.com/lh3/minimap2/archive/refs/tags/v$(minimap2_ver).tar.gz \
		&& tar xzf v$(minimap2_ver).tar.gz \
		&& cd minimap2-$(minimap2_ver) \
		&& make \
		&& mv minimap2 .. \
		&& cd .. \
		&& rm -rf minimap2-* v$(minimap2_ver).tar.gz

mosdepth:
	wget --quiet https://github.com/brentp/mosdepth/releases/download/v$(mosdepth_ver)/mosdepth \
		&& chmod 755 mosdepth

/usr/bin/pandoc:
	apt update \
		&& apt install -y pandoc

readme: data tools
	./create_readme.sh

clean:
	rm -rf genome/chrX.fa github-markdown-toc samtools bwa minimap2 mostdepth tmp.html && eg/clean.sh
