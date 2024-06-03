import os
import subprocess
import numpy as np
import pandas as pd
import multiprocessing
from itertools import repeat
from pyfaidx import Fasta
import pysam
from collections import defaultdict

from seacon.helpers import flatten

# Get chromosome lengths from a reference genome.
def get_chrom_lens(ref_path, chrom_names):
    ref = Fasta(ref_path)
    chrom_lens = {}
    for chrom in chrom_names:
        chrom_lens[chrom] = len(ref[chrom])
    return chrom_lens

# Returns list of tuples containing regions
def get_bins_from_chromlens(chrom_lens, bin_size):
    regions = []
    chrom_names = list(chrom_lens.keys())
    #chrom_names.sort(key=lambda x: int(x[3:]))
    for chrom, tot_len in chrom_lens.items():
        start = 1
        while start < tot_len:
            end = min(start + bin_size - 1, tot_len)
            regions.append((chrom, start, end))
            start += bin_size
    return regions

# Creates bed file containing chrom, start, and end indices of regions from regions list
def write_region_file(out_dir, file_name, regions, stats=None):
    file_path = os.path.join(out_dir, file_name)
    with open(file_path, 'w+') as f:
        for r in regions:
            if stats:
                g,m = stats[(r[0], r[1])]
                f.write(f'{r[0]}\t{r[1]}\t{r[2]}\t{g}\t{m}\n')
            else:
                f.write(f'{r[0]}\t{r[1]}\t{r[2]}\n')
    return file_path

# Read regions from a bed file and returns them as a list of tuples
def read_region_file(region_path):
    regions = []
    with open(region_path) as f:
        for line in f:
            info = line.strip().split('\t')
            chrom, start, end = info[0], info[1], info[2]
            regions.append((chrom, start, end))
    return regions

def read_stats(region_path):
    gc, mapp = [], []
    with open(region_path) as f:
        for line in f:
            info = line.strip().split('\t')
            g,m = float(info[3]), float(info[4])
            gc.append(g)
            mapp.append(m)
    return gc, mapp

# Writes cell names to file
def write_cell_names(path, cell_names):
    with open(path, 'w+') as f:
        for cell in cell_names:
            f.write(f'{cell}\n')

# Reads cell names from file
def read_cell_names(path):
    cell_names = []
    with open(path) as f:
        for line in f:
            cell_names.append(line.strip())
    return cell_names

def get_chrom_names(inp):
    if inp is None:
        chrom_names = [f'chr{i}' for i in range(1, 23)]
    else:
        with open(inp) as f:
            chrom_names = [line.strip() for line in f]
    return chrom_names

        
def get_readcounts_bed(bed_path, bin_ids):
    f = open(bed_path)
    lines = f.readlines()
    f.close()

    counts = {}
    for line in lines:
        info = line.split('\t')
        chrom, start, end, reads = info[0], int(info[1]), int(info[2]), int(info[3])
        counts[bin_ids[chrom][start]] = reads
    return counts

## Unfinished
def get_readcounts_bamfile(bam_path, regions_inp):
    if isinstance(regions_inp, str):
        regions = read_region_file(region_inp)
    else:
        regions = regions_inp
    f = pysam.AlignmentFile(bam_path, 'rb')
    counts = [f.count(region=f'{r[0]}:{r[1]}-{r[2]}') for r in regions]
    return np.array(counts)

def get_readcounts_bamdir(bam_dir, regions, num_processors=1):
    bam_paths = [fname for fname in os.listdir(bam_dir) if os.path.isfile(os.path.join(bam_dir, fname)) and fname[-4:] == '.bam']
    full_bam_paths = [os.path.join(bam_dir, bam_file) for bam_file in bam_paths]
    cell_ids = [name[:-4] for name in bam_paths]

    if num_processors > 1:
        with multiprocessing.Pool(processes=num_processors) as pool:
            counts = pool.starmap(get_readcounts_bamfile, zip(full_bam_paths, repeat(regions)))
    else:
        counts = []
        for bam_file in bam_paths:
            f = pysam.AlignmentFile(os.path.join(bam_dir, bam_file), 'rb')
            cur_counts = [f.count(region=f'{r[0]}:{r[1]}-{r[2]}') for r in regions]
            counts.append(cur_counts)
            f.close()

    readcounts_df = pd.DataFrame(counts, index=cell_ids)
    readcounts_df.index.name = 'cell'
    return readcounts_df, cell_ids

# Collects GC content over regions specified by bed file
def get_gc(ref_path, bed_path, bedtools_path):
    run = subprocess.run([bedtools_path, 'nuc', '-fi', ref_path, '-bed', bed_path], check=True, stdout=subprocess.PIPE, text=True)
    lines = run.stdout.strip().split('\n')
    full_info = [line.split('\t') for line in lines[1:]]
    gc = [float(x[4]) for x in full_info]
    return gc

# Gets mappability info over regions from bed file from bigwig file
def get_mapp(outdir, bed_path, bigwig_path, map_file):
    temp_bed = os.path.join(outdir, 'temp_regions.bed')
    with open(bed_path) as f1, open(temp_bed, 'w+') as f2:
        lines = f1.readlines()
        for i,line in enumerate(lines):
            f2.write(f'{line.strip()}\t{i}\n')

    map_out = os.path.join(outdir, 'mappability.tab')
    run = subprocess.run([bigwig_path, map_file, temp_bed, map_out], check=True, text=True)
    mapp = []
    with open(map_out) as f:
        for line in f:
            score = float(line.strip().split('\t')[-1])
            mapp.append(score)
    os.remove(temp_bed)
    os.remove(map_out)
    return mapp

# Collects data from mappability file text
def read_map(map_path, bin_ids):
    map_scores = {}
    with open(map_path) as f:
        next(f)
        for line in f:
            info = line.strip().split('\t')
            chrom, start, mapp, at, gc, nn = 'chr' + info[0], int(info[1]) - 1, float(info[3]), float(info[4]), float(info[5]), float(info[6])
            map_scores[bin_ids[chrom][start]] = [mapp, at, gc, nn]
    return map_scores

# Filters regions based on gc and mapp
def filter_bins_gc_mapp(regions, gc, mapp, min_gc = 0.2, max_gc = 0.8, min_mapp = 0.9):
    filtered_regions = []
    filtered_stats = {}
    for i,region in enumerate(regions):
        if gc[i] > min_gc and gc[i] < max_gc and mapp[i] > min_mapp:
            filtered_regions.append(region)
            filtered_stats[(region[0], region[1])] = (gc[i], mapp[i])
    return filtered_regions, filtered_stats

def get_chisel_barcodes(filepath, inverse=False):
    df = pd.read_csv(filepath, sep='\t')
    df = df.drop(['REPETITION', 'FILE'], axis=1, errors='ignore')
    df['#CELL'] = df['#CELL'].apply(lambda x: x[:-4])
    df = df.set_index('BARCODE')
    df = df.to_dict()['#CELL']
    if inverse:
        df = {v: k for k, v in df.items()}
    return df

def collect_other_BAF(filepath, full_cell_names):
    df = pd.read_csv(filepath, sep='\t')
    df = df.set_index(['cell', 'chrom', 'start', 'end'])

    alt_cell_names = list(df.index.get_level_values(0).drop_duplicates())
    cell_df = df.xs(alt_cell_names[0])
    regions = [(x[0], x[1]+1, x[2]) for x in cell_df.index]
    df = df.droplevel('end')
    df['BAF'] = df['BAF'].apply(lambda x: min(x, 1-x))
    cell_names = [c for c in full_cell_names if c in alt_cell_names]
    df = flatten(df, cell_names)
    return df, regions

# Read readcounts from file in simple format (header is 'cell' + 0...n, where n is number of bins, 1 row for each cell)
def read_data_flat(path):
    df = pd.read_csv(path, delimiter='\t', index_col=0)
    df.columns = df.columns.map(int)
    return df


# Read readcounts or bafs from file in complete format
def read_data_full(path):
    df = pd.read_csv(path, delimiter='\t')
    df = df.drop(['end'], axis=1)
    df = df.set_index(['cell', 'chrom', 'start'])
    return df