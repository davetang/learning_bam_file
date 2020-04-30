## README

It is often useful to know the start and end coordinates of assembled chromosomes and contigs. We can obtain this information from the MySQL server hosted by the UCSC Genome Browser team.

```bash
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo" > hg19_info.tsv
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg38.chromInfo" > hg38_info.tsv
```

If you set your MySQL config file (`~/.my.cnf`) as

```
[clientucsc]
user=genome
password=
host=genome-mysql.cse.ucsc.edu
```

you can run the following instead:

```bash
mysql --defaults-group-suffix=ucsc -A -e "select chrom, size from hg19.chromInfo" > hg19_info.tsv
mysql --defaults-group-suffix=ucsc -A -e "select chrom, size from hg38.chromInfo" > hg38_info.tsv
```

