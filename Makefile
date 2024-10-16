# Usage:
# make install SRC=/path/to/your/directory
# make liftover VAR_TAB_FILE=test/test.tsv OUT_PREFIX=test_prefix OUT_DIR=test/ SRC=/path/to/your/directory

.PHONY: install help liftover check_vars create_out_dir convert_tsv_to_vcf compress_and_index liftover_variants normalize_variants query_variants cleanup

help:
	@printf "\nINSTALLATION\n\tUsage: make install SRC=/path/to/your/directory\n\tSRC is optional. If SRC is not provided, the files are downloaded to the current directory.\n"
	@printf "\nLIFTOVER\n\tUsage: make liftover VAR_TAB_FILE=test/test.tsv OUT_PREFIX=test_prefix OUT_DIR=test/ SRC=/path/to/your/directory\n\tOUT_DIR is optional. If OUT_DIR is not provided, the output will be stored in the current directory.\n\tSRC is optional. If you used a different directory for installation, then specify the path SRC=/path/to/your/directory\n"

SRC ?= $(PWD)
SRC_DIR := $(shell realpath $(SRC))

install: $(SRC_DIR)/genomes/hg19.fa $(SRC_DIR)/genomes/hg38.fa $(SRC_DIR)/genomes/hg19ToHg38.over.chain.gz $(SRC_DIR)/tools/libexec/bcftools/liftover.so
	@printf "\nUse this genome dir in the liftover shell script:\n$(SRC_DIR)/genomes/\n\n"
	@printf "Use this bcftools in the liftover shell script:\n$(SRC_DIR)/tools/bin/bcftools\n\n"

$(SRC_DIR)/genomes/hg19.fa:
	@printf "Downloading hg19.fa\n"
	@mkdir -p $(SRC_DIR)/genomes
	@wget --no-check-certificate --no-verbose https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz -O $(SRC_DIR)/genomes/hg19.fa.gz
	@gunzip $(SRC_DIR)/genomes/hg19.fa.gz

$(SRC_DIR)/genomes/hg38.fa:
	@printf "Downloading hg38.fa\n"
	@mkdir -p $(SRC_DIR)/genomes
	@wget --no-check-certificate --no-verbose https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -O $(SRC_DIR)/genomes/hg38.fa.gz
	@gunzip $(SRC_DIR)/genomes/hg38.fa.gz

$(SRC_DIR)/genomes/hg19ToHg38.over.chain.gz:
	@printf "Downloading hg19ToHg38.over.chain\n"
	@mkdir -p $(SRC_DIR)/genomes
	@wget --no-check-certificate --no-verbose https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz -O $(SRC_DIR)/genomes/hg19ToHg38.over.chain.gz

$(SRC_DIR)/tools/libexec/bcftools/liftover.so:
	@printf "Installing liftover\n"
	@mkdir -p $(SRC_DIR)/tools
	@wget --no-check-certificate --no-verbose https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 -O $(SRC_DIR)/tools/bcftools-1.19.tar.bz2
	@tar -xjf $(SRC_DIR)/tools/bcftools-1.19.tar.bz2 -C $(SRC_DIR)/tools
	@cd $(SRC_DIR)/tools/bcftools-1.19 && ./configure --prefix=$(SRC_DIR)/tools && make && make install

	@wget --no-check-certificate --no-verbose https://software.broadinstitute.org/software/score/score_1.20-20240927.zip -O $(SRC_DIR)/tools/score_1.20-20240927.zip
	@unzip $(SRC_DIR)/tools/score_1.20-20240927.zip -d $(SRC_DIR)/tools/libexec/bcftools
	
	@wget --no-check-certificate --no-verbose https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2 -O $(SRC_DIR)/tools/htslib-1.19.1.tar.bz2
	@tar -xjf $(SRC_DIR)/tools/htslib-1.19.1.tar.bz2 -C $(SRC_DIR)/tools
	@cd $(SRC_DIR)/tools/htslib-1.19.1 && ./configure --prefix=$(SRC_DIR)/tools && make && make install

	@rm -rf $(SRC_DIR)/tools/bcftools-1.19.tar.bz2 $(SRC_DIR)/tools/score_1.20-20240927.zip $(SRC_DIR)/tools/bcftools-1.19 $(SRC_DIR)/tools/htslib-1.19.1.tar.bz2 $(SRC_DIR)/tools/htslib-1.19.1

SCRIPT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
REFERENCE := $(SRC_DIR)/genomes
TOOLS_DIR := $(SRC_DIR)/tools/bin
PYTHON3_PATH := $(shell command -v python3)

# Check if Python 3 is installed
ifeq ($(PYTHON3_PATH),)
$(error Python 3 is not installed. Please install it and try again.)
endif

# Output directory
OUT_DIR ?= $(OUT_PREFIX)_liftover_results

liftover: check_vars create_out_dir convert_tsv_to_vcf compress_and_index liftover_variants normalize_variants query_variants cleanup

check_vars:
ifndef VAR_TAB_FILE
	$(warning Usage: make liftover VAR_TAB_FILE=test/test.tsv OUT_PREFIX=test_prefix OUT_DIR=test/)
	$(error VAR_TAB_FILE is not set. Please set it to the input tab file.)
endif
ifndef OUT_PREFIX
	$(warning Usage: make liftover VAR_TAB_FILE=test/test.tsv OUT_PREFIX=test_prefix OUT_DIR=test/)
	$(error OUT_PREFIX is not set. Please set it to the output prefix.)
endif

create_out_dir:
	@mkdir -p $(OUT_DIR)

convert_tsv_to_vcf:
	@$(PYTHON3_PATH) $(SCRIPT_DIR)/scripts/convert_tsv_to_vcf.py --tab $(VAR_TAB_FILE) --vcf $(OUT_DIR)/$(OUT_PREFIX).vcf

compress_and_index:
	@$(TOOLS_DIR)/bgzip $(OUT_DIR)/$(OUT_PREFIX).vcf
	@$(TOOLS_DIR)/tabix -p vcf $(OUT_DIR)/$(OUT_PREFIX).vcf.gz

liftover_variants:
	@$(TOOLS_DIR)/bcftools +liftover -Ov -o $(OUT_DIR)/$(OUT_PREFIX)_lo_hg19Tohg38.vcf $(OUT_DIR)/$(OUT_PREFIX).vcf.gz -- -s $(REFERENCE)/hg19.fa -f $(REFERENCE)/hg38.fa -c $(REFERENCE)/hg19ToHg38.over.chain.gz --reject $(OUT_DIR)/$(OUT_PREFIX).rejected.vcf -Ov --write-src

normalize_variants:
	@$(TOOLS_DIR)/bcftools norm --rm-dup exact $(OUT_DIR)/$(OUT_PREFIX)_lo_hg19Tohg38.vcf -Ov -o $(OUT_DIR)/$(OUT_PREFIX)_lo_hg19Tohg38_norm.vcf

query_variants:
	@$(TOOLS_DIR)/bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\tSRC_CHROM=%SRC_CHROM;SRC_POS=%SRC_POS;SRC_ID=%ID;SRC_REF_ALT=%SRC_REF_ALT" $(OUT_DIR)/$(OUT_PREFIX)_lo_hg19Tohg38_norm.vcf > $(OUT_DIR)/$(OUT_PREFIX)_lo_variants.tsv
	@echo "Liftover variants list: $(OUT_DIR)/$(OUT_PREFIX)_lo_variants.tsv"

cleanup:
ifneq ($(CLEANUP),)
ifneq ($(CLEANUP), clean)
	@rm $(OUT_DIR)/$(OUT_PREFIX).vcf* $(OUT_DIR)/$(OUT_PREFIX)_lo_hg19Tohg38.vcf
endif
endif