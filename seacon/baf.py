import os
import subprocess
import multiprocessing
from itertools import repeat
import heapq
import numpy as np
import pandas as pd
import pysam
import vcf
import pickle
import logging

import seacon.reader as reader
from seacon.helpers import index_bins

def create_barcodes(cell_names, temp_chis_dir):
    characters = ['A', 'T', 'C', 'G']
    barcodes = []
    cell_map = {}
    with open(os.path.join(temp_chis_dir, 'barcodes.info.tsv'), 'w+') as f:
        headers = ['#CELL', 'BARCODE']
        f.write('\t'.join(headers) + '\n')
        for cell in cell_names:
            cur_bc = ''.join(np.random.choice(characters, size=12))
            while cur_bc in barcodes:
                cur_bc = np.random.choice(characters, size=12)
            barcodes.append(cur_bc)
            cell_map[cell] = cur_bc
            f.write(f'{cell}.bam\t{cur_bc}\n')
    return cell_map

def gen_totalcount_inp(cell_names, cell_map, readcounts, temp_chis_dir):
    with open(os.path.join(temp_chis_dir, 'rdr', 'total.tsv'), 'w+') as f:
        for cell in cell_names:
            f.write(f'{cell_map[cell]}\t{str(sum(readcounts.loc[cell]))}\n')

def gen_rdr_inp(cell_names, cell_map, bin_coords, norm_counts, readcounts, RDR, temp_chis_dir):
    with open(os.path.join(temp_chis_dir, 'rdr', 'rdr.tsv'), 'w+') as f:
        for i,b in enumerate(bin_coords):
            chrom, start, end = b[0], b[1]-1, b[2]
            for cell in cell_names:
                f.write(f'{chrom}\t{start}\t{end}\t{cell_map[cell]}\t{norm_counts[i]}\t{readcounts.loc[cell][i]}\t{RDR.loc[cell][i]}\n')

def get_phased_counts_cell(cell, chrom_names, vcf_path, bam_dir, temp_chis_dir):
    logging.info(f'Collecting phased counts for cell {cell}')
    bam_path = os.path.join(bam_dir, f'{cell}.bam')
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    vcf_reader = vcf.Reader(filename=vcf_path)

    counts = {chrom: [] for chrom in chrom_names}
    for r in vcf_reader:
        chrom = r.CHROM
        pos = r.POS
        ref_base = r.REF
        alt_bases = r.ALT

        ref_count = 0
        alt_count = 0

        for read in bamfile.fetch(chrom, pos - 1, pos):
            positions = read.get_reference_positions(full_length=True)
            if pos - 1 in positions:
                read_pos = positions.index(pos - 1)
                if read_pos is not None and read_pos < len(read.query_sequence):
                    base = read.query_sequence[read_pos]
                    if base == ref_base:
                        ref_count += 1
                    elif base in alt_bases:
                        alt_count += 1
        
        if ref_count + alt_count > 0:
            counts[chrom].append((pos, ref_count, alt_count, cell))

    with open(os.path.join(temp_chis_dir, 'baf', f'{cell}_temp.pkl'), 'wb') as f:
        pickle.dump(counts, f)

def merge_phased_counts(cell_names, chrom_names, cell_map, temp_chis_dir):
    sample_counts = []
    for cell in cell_names:
        with open(os.path.join(temp_chis_dir, 'baf', f'{cell}_temp.pkl'), 'rb') as f:
            cur_counts = pickle.load(f)
            sample_counts.append(cur_counts)
    with open(os.path.join(temp_chis_dir, 'baf', 'baf.tsv'), 'w+') as f:
        for chrom in chrom_names:
            chrom_counts = [x[chrom] for x in sample_counts]
            h = heapq.merge(*chrom_counts, key = lambda x: x[0])
            for ent in h:
                pos, A, B, cell = ent
                f.write(f'{chrom}\t{pos}\t{cell_map[cell]}\t{A}\t{B}\n')
    for cell in cell_names:
        os.remove(os.path.join(temp_chis_dir, 'baf', f'{cell}_temp.pkl'))
    
def gen_baf_inp(cell_names, chrom_names, cell_map, vcf_path, bam_dir, temp_chis_dir, num_processors):
    bam_paths = [os.path.join(bam_dir, f'{cell}.bam') for cell in cell_names]
    if num_processors > 1:
        with multiprocessing.Pool(processes=num_processors) as pool:
            pool.starmap(get_phased_counts_cell, zip(cell_names, repeat(chrom_names), repeat(vcf_path), repeat(bam_dir), repeat(temp_chis_dir)))
            
    else:
        for cell in cell_names:
            get_phased_counts_cell(cell, chrom_names, vcf_path, bam_dir, temp_chis_dir)
    merge_phased_counts(cell_names, chrom_names, cell_map, temp_chis_dir)

def call_chisel_combo(chisel_src, temp_chis_dir, block_size, num_processors):
    tot_path = os.path.join(temp_chis_dir, 'rdr', 'total.tsv')
    rdr_path = os.path.join(temp_chis_dir, 'rdr', 'rdr.tsv')
    baf_path = os.path.join(temp_chis_dir, 'baf', 'baf.tsv')
    out_path = os.path.join(temp_chis_dir, 'combo', 'combo.tsv')
    cmd = ["conda", "run", "-n", "SEACON_chisel", "python2.7", os.path.join(chisel_src, "Combiner.py"), "-r", rdr_path, "-b", baf_path, "-j", str(num_processors), "-k", f'{block_size}kb', "-l", tot_path]
    with open(out_path, 'w+') as f:
        subprocess.run(cmd, stdout=f, text=True)

def call_chisel_caller(chisel_src, temp_chis_dir, max_CN, num_processors):
    combo_path = os.path.join(temp_chis_dir, 'combo', 'combo.tsv')
    out_path = os.path.join(temp_chis_dir, 'calls', 'calls.tsv')
    cmd = ["conda", "run", "-n", "SEACON_chisel", "python2.7", os.path.join(chisel_src, "Caller.py"), combo_path, "-P", str(max_CN), "-j", str(num_processors)]
    with open(out_path, 'w+') as f:
        subprocess.run(cmd, stdout=f, text=True)

def chisel_baf_helper(cell_names, chrom_names, bin_coords, norm_counts, readcounts_df, RDR_df, bam_dir, vcf_path, out_dir, temp_chis_dir, num_processors):
    if os.path.exists(os.path.join(temp_chis_dir, 'barcodes.info.tsv')):
        cell_map = reader.get_chisel_barcodes(os.path.join(temp_chis_dir, 'barcodes.info.tsv'), inverse=True)
    else:
        cell_map = create_barcodes(cell_names, temp_chis_dir)
    
    new_dirs = ['baf', 'calls', 'combo', 'plots', 'rdr']
    for d in new_dirs:
        if not os.path.isdir(os.path.join(temp_chis_dir, d)):
            os.mkdir(os.path.join(temp_chis_dir, d))

    gen_totalcount_inp(cell_names, cell_map, readcounts_df, temp_chis_dir)
    gen_rdr_inp(cell_names, cell_map, bin_coords, norm_counts, readcounts_df, RDR_df, temp_chis_dir)
    gen_baf_inp(cell_names, chrom_names, cell_map, vcf_path, bam_dir, temp_chis_dir, num_processors)

def convert_baf(out_dir, temp_chis_dir, cell_names, bin_coords):
    cell_map = reader.get_chisel_barcodes(os.path.join(temp_chis_dir, 'barcodes.info.tsv'))

    df = pd.read_csv(os.path.join(temp_chis_dir, 'combo', 'combo.tsv'), sep='\t', names=['chrom', 'start', 'end', 'cell', 'normcount', 'readcount', 'RDR', 'Acount', 'Bcount', 'BAF'])
    df = df.drop(['end', 'normcount', 'readcount', 'RDR', 'Acount', 'Bcount'], axis=1)
    df['BAF'] = df['BAF'].apply(lambda x: min(x, 1-x))
    df['start'] = df['start'] + 1
    df['cell'] = df['cell'].apply(lambda x: cell_map[x])
    df['id'] = df['chrom'] + '-' + df['start'].astype(str)
    df = df.pivot(index='cell', columns='id', values='BAF')
    order = [f'{x[0]}-{x[1]}' for x in bin_coords]
    df = df[order]
    df = df.rename(columns = {x: i for i,x in enumerate(order)})
    df = df.reindex(cell_names)
    df = df.round(5)
    df.to_csv(os.path.join(out_dir, 'BAF.tsv'), sep='\t')