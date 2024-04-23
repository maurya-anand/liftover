#!/bin/bash

# liftover.hg19.to.hg38.sh
#
# This script converts of genomic coordinates from the Human Genome version 19 (hg19) to version 38 (hg38) using bcftool's liftover plugin.
# It takes as input a tab-separated values (TSV) file containing variant information with columns specifying chromosome (chr), position (pos), identifier (id), reference allele (ref), and alternate allele (alt). 
#
# Usage:
#   bash liftover.hg19.to.hg38.sh <input_tab_file> <output_prefix>
#
# Parameters:
#   input_tab_file: The input file which should be a tab-delimited file with the following columns: chr, pos, id, ref, and alt.
#   output_prefix: The prefix for the output file. The script will add a .vcf extension to this prefix to create the output file name.
#
# Example:
#   bash liftover.hg19.to.hg38.sh test.tsv test_set
#
# This will create the output file test_set_lo_variants.tsv in the test_set_liftover_results directory.

VAR_TAB_FILE=$1
OUT_PREFIX=$2
OUT_DIR=$3
CLEANUP=$4

SCRIPT_DIR=$(dirname "$0")
REFERENCE="genomes"
TOOLS_DIR=${SCRIPT_DIR}/tools/bin/

# Check if Python 3 is installed
if command -v python3 &> /dev/null
then
    PYTHON3_PATH=$(command -v python3)
else
    echo "Python 3 is not installed. Please install it and try again."
    exit 1
fi

# If OUT_DIR is provided, create a subdirectory under it
if [ -n "$OUT_DIR" ]; then
    OUT_DIR=${OUT_DIR}/${OUT_PREFIX}_liftover_results
else
    OUT_DIR=${OUT_PREFIX}_liftover_results
fi

mkdir -p ${OUT_DIR}

# convert tsv to vcf
$PYTHON3_PATH ${SCRIPT_DIR}/scripts/convert_tsv_to_vcf.py --tab ${VAR_TAB_FILE} --vcf ${OUT_DIR}/${OUT_PREFIX}.vcf

# compress and generate vcf index
${TOOLS_DIR}/bgzip ${OUT_DIR}/${OUT_PREFIX}.vcf
${TOOLS_DIR}/tabix -p vcf ${OUT_DIR}/${OUT_PREFIX}.vcf.gz

# liftover variants
${TOOLS_DIR}/bcftools +liftover -Ov -o ${OUT_DIR}/${OUT_PREFIX}_lo_hg19Tohg38.vcf ${OUT_DIR}/${OUT_PREFIX}.vcf.gz -- -s ${REFERENCE}/hg19.fa -f ${REFERENCE}/hg38.fa -c ${REFERENCE}/hg19ToHg38.over.chain.gz --reject ${OUT_DIR}/${OUT_PREFIX}.rejected.vcf -Ov --write-src

${TOOLS_DIR}/bcftools norm --rm-dup exact ${OUT_DIR}/${OUT_PREFIX}_lo_hg19Tohg38.vcf -Ov -o ${OUT_DIR}/${OUT_PREFIX}_lo_hg19Tohg38_norm.vcf

${TOOLS_DIR}/bcftools query -f "%CHROM\t%POS\tsrc:%SRC_CHROM,%SRC_POS,%ID,%SRC_REF_ALT\t%REF\t%ALT" ${OUT_DIR}/${OUT_PREFIX}_lo_hg19Tohg38_norm.vcf | sed 's/,/:/g' > ${OUT_DIR}/${OUT_PREFIX}_lo_variants.tsv

echo "Liftover variants list: ${OUT_DIR}/${OUT_PREFIX}_lo_variants.tsv"

if [ "${CLEANUP}" = "clean" ]; then
    rm ${OUT_DIR}/${OUT_PREFIX}.vcf* ${OUT_DIR}/${OUT_PREFIX}_lo_hg19Tohg38.vcf
fi