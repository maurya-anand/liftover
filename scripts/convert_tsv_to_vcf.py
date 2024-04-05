import os
import argparse
import pandas as pd

"""
convert_to_vcf.py

This script is used to convert a tab-delimited file with specific columns into a VCF (Variant Call Format) file.

Functions:
----------
create_vcf(input_file: str, output_file: str) -> None:
    This function takes an input file and an output file as arguments. The input file should be a tab-delimited file with the following columns: chr, pos, id, ref, and alt. The function reads the input file, processes the data, and writes it to the output file in VCF format. If any of the required data (chr, pos, ref, alt) is missing in a row, it prints an error message and writes a VCF line with empty fields.

Main Execution:
---------------
The script takes two command-line arguments:
    --tab: The input file which contains the following columns: chr, pos, id, ref, and alt. This argument is required.
    --vcf: The name of the output VCF file (including the file path). This argument is required.

Example usage:
    python3 convert_tsv_to_vcf.py --tab temp_coordinates.txt --vcf temp_coordinates.vcf
"""

def create_vcf(input_file, output_file):
    """
    Convert a tab-delimited file into a VCF (Variant Call Format) file.

    Parameters:
    input_file (str): The input file which should be a tab-delimited file with the following columns: chr, pos, id, ref, and alt.
    output_file (str): The output file where the VCF data will be written.

    Returns:
    None
    """
    sample_name = os.path.splitext(os.path.basename(output_file))[0]

    # Read the input file into a DataFrame
    df = pd.read_csv(input_file, sep='\t', names=['chr', 'pos', 'id', 'ref', 'alt'])

    # Replace missing ids with "."
    df['id'] = df['id'].fillna('.')

    # Sort the DataFrame by chromosome and position
    df.sort_values(['chr', 'pos'], inplace=True) 

    # Write the DataFrame to the output file in VCF format
    with open(output_file, 'w') as vcf_out:
        vcf_out.write("##fileformat=VCFv4.2\n")
        vcf_out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")
        for _, row in df.iterrows():
            vcf_out.write(f"{row['chr']}\t{row['pos']}\t{row['id']}\t{row['ref']}\t{row['alt']}\t.\t.\t.\t.\t.\n")
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Convert a tab-delimited file into a VCF (Variant Call Format) file.',
        epilog='Ensure that the input file contains the following columns: chr, pos, id, ref, and alt.\n\nExample usage: python3 convert_tsv_to_vcf.py --tab temp_coordinates.txt --vcf temp_coordinates.vcf',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--tab', help='The input file which should be a tab-delimited file with the following columns: chr, pos, id, ref, and alt.', required=True)
    parser.add_argument('--vcf', help='The output file where the VCF data will be written.', required=True)
    args = parser.parse_args()
    create_vcf(args.tab, args.vcf)
