all: data tools readme
data: genome/chrX.fa
tools: github-markdown-toc samtools miniconda3 pandoc miniconda3/bin/R miniconda3/lib/R/library/rmarkdown

genome/chrX.fa:
	bunzip2 -c eg/chrX.fa.bz2 > genome/chrX.fa

github-markdown-toc:
	git clone https://github.com/ekalinin/github-markdown-toc.git

samtools:
	wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 && tar xjf samtools-1.13.tar.bz2 && cd samtools-1.13 && ./configure && make && mv samtools .. && cd .. && rm -rf samtools-1.13.tar.bz2 samtools-1.13

miniconda3:
	wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-Linux-x86_64.sh && bash Miniconda3-py39_4.9.2-Linux-x86_64.sh -b -p miniconda3 && miniconda3/bin/conda update -y -n base -c defaults conda && rm Miniconda3-py39_4.9.2-Linux-x86_64.sh

pandoc: miniconda3
	miniconda3/bin/conda install -y -c conda-forge pandoc

miniconda3/bin/R: miniconda3
	miniconda3/bin/conda install -y -c conda-forge r-base

miniconda3/lib/R/library/rmarkdown: miniconda3/bin/R
	miniconda3/bin/Rscript -e "Sys.setenv(PATH=paste0(Sys.getenv('PATH'), ':', getwd(), '/miniconda3/bin')); install.packages('rmarkdown', repos='http://cran.us.r-project.org')"

readme: genome/chrX.fa github-markdown-toc samtools miniconda3 pandoc miniconda3/bin/R miniconda3/lib/R/library/rmarkdown
	./create_readme.sh

clean:
	rm -rf genome/chrX.fa github-markdown-toc samtools miniconda3 && eg/clean.sh
