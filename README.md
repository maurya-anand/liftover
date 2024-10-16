# liftover

![GitHub Workflow Status (with event)](https://img.shields.io/github/actions/workflow/status/maurya-anand/liftover/publish.yml?style=flat-square)&nbsp;&nbsp;
![GitHub release (with filter)](https://img.shields.io/github/v/release/maurya-anand/liftover?style=flat-square)&nbsp;&nbsp;

This script performs the conversion of genomic coordinates from Human Genome version 19 (hg19) to version 38 (hg38) using the liftover plugin from bcftools.

It takes a tab-separated values (TSV) file as input. This file should contain variant information with columns denoting chromosome (chr), position (pos), identifier (id), reference allele (ref), and alternate allele (alt).

> [!TIP]
> Advantage of using `BCFtools/liftover` over other tools: [Comparison of features and limitations across `BCFtools/liftover`, `Transanno/liftover`, `Genozip/DVCF`, `GenomeWrap`, `Picard/LiftoverVcf` and `CrossMap/VCF`](https://academic.oup.com/view-large/438467641)
>
> Genovese, Giulio, _et al._ "**BCFtools/Liftover: An Accurate and Comprehensive Tool to Convert Genetic Variants across Genome Assemblies.**" _Bioinformatics_, vol. 40, no. 2, 2024, <https://doi.org/10.1093/bioinformatics/btae038>.

## Quickstart

```
cd liftover
make help

INSTALLATION
        Usage: make install SRC=/path/to/your/directory
        SRC is optional. If SRC is not provided, the files are downloaded to the current directory.

LIFTOVER
        Usage: make liftover VAR_TAB_FILE=test/test.tsv OUT_PREFIX=test_prefix OUT_DIR=test/ SRC=/path/to/your/directory
        OUT_DIR is optional. If OUT_DIR is not provided, the output will be stored in the current directory.
        SRC is optional. If you used a different directory for installation, then specify the path SRC=/path/to/your/directory
```

## Installation

> [!IMPORTANT]
> Please note that this script requires `python3` and the `pandas` library. However, if you are using the Docker image, these requirements are already included in the image.

To install the necessary tools, which include `bcftools`, `htslib`, and the `liftover plugin`, and to download the required genome files (hg19.fa, hg38.fa and hg19ToHg38.over.chain.gz), execute the following command:

``` bash
cd liftover
make install
```

This command will also create a directory named `genomes` and store the downloaded genome files in it. If `SRC` is not provided, the files are downloaded to the current directory.

## Liftover

### Using the Makefile (recommended)

``` bash
cd liftover
make liftover VAR_TAB_FILE=/path/to/your/input_tab_file.tsv OUT_PREFIX=output_prefix
```

### Parameters

- `VAR_TAB_FILE`: The input file should be a tab-delimited file containing the following columns: chr, pos, id, ref, and alt.
- `OUT_PREFIX`: Prefix for the output file.
- `OUT_DIR`: Output directory. (optional)
- `SRC`: Path to the directory if `make install SRC=/path/to/your/directory` was used during installation. (optional)

### Using the Bash script

After installing the necessary tools and downloading the required genome files using the steps provided in the "Installation" section, you can run the script as follows:

``` bash
cd liftover
bash liftover_hg19_to_hg38.sh <input_tab_file> <output_prefix> <output_directory>
```

### Parameters

- `input_tab_file`: The input file should be a tab-delimited file containing the following columns: chr, pos, id, ref, and alt.
- `output_prefix`: Prefix for the output file.


### Using Docker

This package is also available as a Docker image, which can be pulled from the GitHub Container Registry using the following command:

``` bash
docker pull ghcr.io/maurya-anand/liftover
```

You can run the Docker image as follows:

``` bash
cd liftover
docker run -v $(pwd):/data -it ghcr.io/maurya-anand/liftover /scripts/liftover_hg19_to_hg38.sh /data/test/test.tsv docker_test /data
```

The output from the above command will be saved in a directory named `docker_test_liftover_results` in your current directory.

In this command, `-v $(pwd):/data` mounts your current directory to the `/data` directory in the Docker container.

Replace `/data/test/test.tsv` with your input file, `docker_test` with your desired output prefix, and `/data` with your desired output directory.

### Example usage

> [!WARNING]
> The script assumes that the liftover plugin and reference genomes are located in specific directories generated via `make install` command. These paths may need to be adjusted based on your specific setup. The script also assumes that the input file is properly formatted and contains the necessary columns.

#### Input file

`cat test/test.tsv`

|     |        |              |     |     |
|-----|--------|--------------|-----|-----|
| 1   | 818046 | .            | T   | C   |
| 2   | 265023 | .            | C   | A   |
| 3   | 361463 | 3:361463:G:T | G   | T   |

#### Execution

``` bash
cd liftover
make liftover VAR_TAB_FILE=test/test.tsv OUT_PREFIX=test_set
```

This command will create a directory named `<prefix>_liftover_results`, where `<prefix>` is the provided prefix, and save the output in a TSV file named `<prefix>_lo_variants.tsv`.

Each variant is represented as a row with the following columns:

- `chr`: The chromosome where the variant is located in the target genome assembly.
- `pos`: The position of the variant on the chromosome in the target genome assembly.
- `ref`: The reference allele in the target genome assembly.
- `alt`: The alternate allele in the target genome assembly.
- `source`: Additional information about the variant in the original genome assembly. This includes:
  - `SRC_CHROM`: The chromosome where the variant is located in the original genome assembly.
  - `SRC_POS`: The position of the variant on the chromosome in the original genome assembly.
  - `SRC_ID`: The identifier of the variant in the original genome assembly.
  - `SRC_REF_ALT`: The reference and alternate alleles in the original genome assembly.

#### Output file

`cat test_set_liftover_results/test_set_lo_variants.tsv`

|      |        |     |     |                                                                |
|-------|---------|------|------|----------------------------------------|
| chr1 | 882666 | T   | C   | SRC_CHROM=1;SRC_POS=818046;SRC_ID=.;SRC_REF_ALT=T,C            |
| chr2 | 265023 | C   | A   | SRC_CHROM=2;SRC_POS=265023;SRC_ID=.;SRC_REF_ALT=C,A            |
| chr3 | 319780 | G   | T   | SRC_CHROM=3;SRC_POS=361463;SRC_ID=3:361463:G:T;SRC_REF_ALT=G,T |

## Steps

1. Convert the tab-delimited file to a VCF file using the convert_tsv_to_vcf.py Python script.
1. Compress the VCF file and generate the index.
1. Perform a liftover operation to convert the variants from one genome build to another using bcftools liftover plugin.
1. Normalize the VCF file using bcftools.
1. Generate a mapping table from the old reference to the new reference using bcftools.

## Components

- **Tools**
  - [htslib v1.19.1](https://www.htslib.org/)
  - [bcftools v1.19](https://doi.org/10.1093/gigascience/giab008)
  - [liftover (plugin) v1.20](https://github.com/freeseek/score)
- **Python3 package**
  - pandas
