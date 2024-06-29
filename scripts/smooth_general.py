import os
import numpy as np
import pandas as pd
import sys
import argparse

sys.path.append("SEACON/seacon")

import seacon
import seacon.reader as reader
from seacon.segment import ensemble_cluster_by_frequency_filter

def get_cbs_bkpts(filepath):
    bkpts = []
    cell_names = []
    with open(filepath) as f:
        for line in f:
            info = line.strip().split('\t')
            cell_names.append(info[0])
            cur_bkpt = [int(x) for x in info[1:-1]]
            bkpts.append(cur_bkpt)
    local_bkpts = dict(zip(cell_names, bkpts))
    return local_bkpts

def get_cn_bkpts(CNs, cells):
    n = CNs.shape[1]
    bkpts = {}
    for cell in cells:
        cell_df = CNs.loc[cell]
        bkpts[cell] = []
        for i in range(1,n):
            if cell_df[i] != cell_df[i - 1]:
                bkpts[cell].append(i)
    return bkpts

def get_upper_bkpts(CNs, bkpts, cells, filt):
    n = CNs.shape[1]
    bkpt_counts = {b: 0 for b in range(n)}
    for b in range(1, n):
        c = 0
        for cell in cells:
            if b in bkpts[cell]:
                bkpt_counts[b] += 1


    new_bkpts = {cell: [] for cell in cells}
    for cell in cells:
        for b in bkpts[cell]:
            if bkpt_counts[b] >= filt:
                new_bkpts[cell].append(b)
    return new_bkpts


def smooth(CNs, bkpts, cells):
    n = CNs.shape[1]
    new_CNs = {cell: [] for cell in cells}
    for cell in cells:
        print(cell)
        cur_CNs = list(CNs.loc[cell])
        cur_bkpts = [0] + list(bkpts[cell]) + [n]
        for i in range(len(cur_bkpts) - 1):
            s,e = cur_bkpts[i], cur_bkpts[i+1]
            l = e - s
            seg = cur_CNs[s:e]
            c = max(seg,key=seg.count)
            for j in range(l):
                new_CNs[cell].append(c)
    return new_CNs

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--cn-path', type=str, required=True, help='Path to input CNPs.')
    parser.add_argument('-b', '--bkpt-path', type=str, required=True, help='Path to input CBS breakpoints.')
    parser.add_argument('-o', '--out-dir', type=str, default='out1', help='Output directory.')
    parser.add_argument('--upper-filter', type=int, default=5, help='Upper bound for breakpoint filter.')
    arguments = parser.parse_args()
    return arguments

def main(cn_path, bkpt_path, out_dir, filt):
    args = parse_args()
    cn_path = args.cn_path
    bkpt_path = args.bkpt_path
    out_dir = args.out_dir
    filt = args.upper_filter

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    CNs = reader.read_data_flat(cn_path)
    n = CNs.shape[1]
    local_bkpts = get_cbs_bkpts(bkpt_path)
    cells = list(set(CNs.index).intersection(set(local_bkpts.keys())))
    cn_bkpts = get_cn_bkpts(CNs, cells)
    upper_bkpts = get_upper_bkpts(CNs, cn_bkpts, cells, filt)
    ensemble_bkpts = ensemble_cluster_by_frequency_filter(local_bkpts, cn_bkpts, upper_bkpts)
    new_CNs = smooth(CNs, ensemble_bkpts, cells)
    with open(os.path.join(out_dir, 'flat_profiles.tsv'), 'w+') as f:
        f.write('\t'.join([str(x) for x in range(n)]) + '\n')
        for cell in cells:
            f.write(cell + '\t' + '\t'.join([str(x) for x in new_CNs[cell]]) + '\n')

if __name__ == '__main__':
    main()