# liftover.hg19.to.hg38.sh

This script is used to convert a tab-delimited file into a VCF (Variant Call Format) file, compress it, generate an index, and perform a liftover operation to convert the variants from one genome build to another.

# Dependencies

The script requires `python3` with the `pandas` library installed. 

To install the necessary tools, including `bcftools`, `htslib` and the `liftover plugin`, and to download the required genome files (hg19.fa, hg38.fa and hg19ToHg38.over.chain.gz), run the following command:

```bash
make install
```
This command will also create a directory named 'genomes' and store the downloaded genome files in it.

# Usage

```bash
bash liftover.hg19.to.hg38.sh <input_tab_file> <output_prefix>
```

## Parameters

- `input_tab_file`: The input file which should be a tab-delimited file with the following columns: chr, pos, id, ref, and alt.
- `output_prefix`: The prefix for the output file. The script will add a .vcf extension to this prefix to create the output file name.

## Example

```bash
bash liftover.hg19.to.hg38.sh test.tsv test_out
```

This will create a VCF file named test_out.vcf from the input file test.tsv, compress it, generate an index, and perform a liftover operation to create a new VCF file named temp_vcf_lo_hg19Tohg38.vcf.

# Steps
1. Convert the tab-delimited file to a VCF file using the convert_tsv_to_vcf.py Python script.
1. Compress the VCF file and generate the index.
1. Perform a liftover operation to convert the variants from one genome build to another using bcftools liftover plugin.
1. Normalize the VCF file using bcftools.
1. Generate a mapping table from the old reference to the new reference using bcftools.

# Components
- **Tools**
  - bgzip
  - [tabix](https://doi.org/10.1093/bioinformatics/btq671)
  - [bcftools](https://doi.org/10.1093/gigascience/giab008)
  - [liftover (plugin)](https://github.com/freeseek/score)
- **Python3 packages**
  - pandas

# Notes
- The script assumes that the liftover plugin and reference genomes are located in specific directories. These paths may need to be adjusted based on your specific setup.
- The script also assumes that the input file is properly formatted and contains the necessary columns.