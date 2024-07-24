"""
This is a single utility script to convert an excel file to a tsv file.
"""
import pandas as pd

def excel2tsv(excel_file, sheet_name, tsv_file):
    df = pd.read_excel(excel_file, sheet_name=sheet_name)
    df.to_csv(tsv_file, sep='\t', index=False)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Convert an excel file to a tsv file.')
    parser.add_argument('--excel_file', '-e', help='The input excel file.')
    parser.add_argument('--sheet_name', '-s', help='The sheet name in the excel file.')
    parser.add_argument('--tsv_file', '-t', help='The output tsv file.')
    args = parser.parse_args()
    excel2tsv(args.excel_file, args.sheet_name, args.tsv_file)