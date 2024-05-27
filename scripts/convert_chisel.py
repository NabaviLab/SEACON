import pandas as pd
import numpy as np
import argparse
import os

# Intended for barcodes.info.tsv
def get_chisel_barcodes(filepath):
    df = pd.read_csv(filepath, sep='\t')
    df = df.drop(['REPETITION', 'FILE'], axis=1)
    df['#CELL'] = df['#CELL'].apply(lambda x: x[:-4])
    df = df.set_index('BARCODE')
    df = df.to_dict()['#CELL']
    return df

# Read existing BAF file from chisel, intended for combo.tsv
def collect_chisel_BAF(filepath, cellmap = None):
    df = pd.read_csv(filepath, sep='\t', names=['chrom', 'start', 'end', 'barcode', 'normcount', 'readcount', 'RDR', 'Acount', 'Bcount', 'BAF'])
    df = df.drop(['normcount', 'readcount', 'RDR', 'Acount', 'Bcount'], axis=1)
    df = df.set_index(['chrom', 'start', 'end', 'barcode'])
    if cellmap:
        df = df.rename(index=cellmap)
    df.index.names = ['chrom', 'start', 'end', 'cell']
    df['BAF'] = df['BAF'].apply(lambda x: min(x, 1-x))
    return df

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--in-path', type=str, default=None, help='Path to directory created by CHISEL.')
parser.add_argument('-o', '--out-path', type=str, default=None, help='Path to converted BAF file.')
args = parser.parse_args()


cm = get_chisel_barcodes(os.path.join(args.in_path, 'barcodes.info.tsv'))
df = collect_chisel_BAF(os.path.join(args.in_path, 'combo', 'combo.tsv'), cellmap = cm)
df.to_csv(args.out_path, sep='\t')