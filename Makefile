# This Makefile downloads and extracts genome files and a bcftools plugin.
# The files are downloaded to a directory specified by the `dir` variable.
# If `dir` is not provided, the files are downloaded to the current directory.

# Usage:
# make install dir=/path/to/your/directory

.PHONY: install help

dir ?= $(PWD)

install: $(dir)/genomes/hg19.fa $(dir)/genomes/hg38.fa $(dir)/genomes/hg19ToHg38.over.chain.gz $(dir)/tools/libexec/bcftools/liftover.so
	@printf "\nUse this genome dir in the liftover shell script:\n$(dir)/genomes/\n\n"
	@printf "Use this bcftools in the liftover shell script:\n$(dir)/tools/bin/bcftools\n\n"

$(dir)/genomes/hg19.fa:
	@mkdir -p $(dir)/genomes
	@wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz -O $(dir)/genomes/hg19.fa.gz
	@gunzip $(dir)/genomes/hg19.fa.gz

$(dir)/genomes/hg38.fa:
	@mkdir -p $(dir)/genomes
	@wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -O $(dir)/genomes/hg38.fa.gz
	@gunzip $(dir)/genomes/hg38.fa.gz

$(dir)/genomes/hg19ToHg38.over.chain.gz:
	@mkdir -p $(dir)/genomes
	@wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz -O $(dir)/genomes/hg19ToHg38.over.chain.gz

$(dir)/tools/libexec/bcftools/liftover.so:
	@mkdir -p $(dir)/tools
	@wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 -O $(dir)/tools/bcftools-1.19.tar.bz2
	@tar -xjf $(dir)/tools/bcftools-1.19.tar.bz2 -C $(dir)/tools
	@cd $(dir)/tools/bcftools-1.19 && ./configure --prefix=$(dir)/tools && make && make install

	@wget https://software.broadinstitute.org/software/score/score_1.19-dev.zip -O $(dir)/tools/score_1.19-dev.zip
	@unzip $(dir)/tools/score_1.19-dev.zip -d $(dir)/tools/libexec/bcftools
	
	@wget https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2 -O $(dir)/tools/htslib-1.19.1.tar.bz2
	@tar -xjf $(dir)/tools/htslib-1.19.1.tar.bz2 -C $(dir)/tools
	@cd $(dir)/tools/htslib-1.19.1 && ./configure --prefix=$(dir)/tools && make && make install

	@rm -rf $(dir)/tools/bcftools-1.19.tar.bz2 $(dir)/tools/score_1.19-dev.zip $(dir)/tools/bcftools-1.19 $(dir)/tools/htslib-1.19.1.tar.bz2 $(dir)/tools/htslib-1.19.1

help:
	@echo "Usage: make install dir=/path/to/your/directory"
	@echo "If `dir` is not provided, the files are downloaded to the current directory."