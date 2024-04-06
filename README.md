# liftover.hg19.to.hg38.sh

This script converts of genomic coordinates from the Human Genome version 19 (hg19) to version 38 (hg38) using bcftool's liftover plugin. It takes as input a tab-separated values (TSV) file containing variant information with columns specifying chromosome (chr), position (pos), identifier (id), reference allele (ref), and alternate allele (alt).

## Installation

> [!IMPORTANT]  
> The script requires `python3` with the `pandas` library installed.

To install the necessary tools, including `bcftools`, `htslib` and the `liftover plugin`, and to download the required genome files (hg19.fa, hg38.fa and hg19ToHg38.over.chain.gz), run the following command:

```bash
cd liftover
make install
```

This command will also create a directory named `genomes` and store the downloaded genome files in it.

## Usage

```bash
cd liftover
bash liftover.hg19.to.hg38.sh <input_tab_file> <output_prefix>
```

### Parameters

- `input_tab_file`: The input file should be a tab-delimited file containing the following columns: chr, pos, id, ref, and alt.
- `output_prefix`: Prefix for the output file.

### Example usage

#### Input file

```bash
cat test/test.tsv
1	818046	.	T	C
2	265023	.	C	A
3	361463	.	G	T
```

#### Execution

```bash
cd liftover
bash liftover.hg19.to.hg38.sh test/test.tsv test_set
```

This command will create a directory named `<prefix>_liftover_results`, where <prefix> is the provided prefix, and save the output in a TSV file named `<prefix>_lo_variants.tsv`.
The TSV file will contain the updated genome coordinates in the following order: chr, pos, id, ref, and alt.
The id column will retain the original genomic coordinates from the input file.

#### Output file

```bash
cat test_set_liftover_results/test_set_lo_variants.tsv
chr1	882666	src:1:818046:T:C	T	C
chr2	265023	src:2:265023:C:A	C	A
chr3	319780	src:3:361463:G:T	G	T

```

## Steps

1. Convert the tab-delimited file to a VCF file using the convert_tsv_to_vcf.py Python script.
1. Compress the VCF file and generate the index.
1. Perform a liftover operation to convert the variants from one genome build to another using bcftools liftover plugin.
1. Normalize the VCF file using bcftools.
1. Generate a mapping table from the old reference to the new reference using bcftools.

## Components

- **Tools**
  - bgzip
  - [tabix](https://doi.org/10.1093/bioinformatics/btq671)
  - [bcftools](https://doi.org/10.1093/gigascience/giab008)
  - [liftover (plugin)](https://github.com/freeseek/score)
- **Python3 packages**
  - pandas

## Notes

- The script assumes that the liftover plugin and reference genomes are located in specific directories. These paths may need to be adjusted based on your specific setup.
- The script also assumes that the input file is properly formatted and contains the necessary columns.
