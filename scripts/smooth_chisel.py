import os
import numpy as np
import pandas as pd
import argparse

import seacon
import seacon.reader as reader
from seacon.helpers import index_bins, get_bins_per_chrom
from seacon.segment import parse_local_bkpts, get_bkpts, filter_breakpoint, ensemble_cluster_by_frequency_filter, get_best_pair
from seacon.caller import get_segs, get_cluster_CNs, allele_caller

def get_info(loc_dir, bin_coords, cell_names):
    nbins = len(bin_coords)
    order = [f'{x[0]}-{x[1]}' for x in bin_coords]
    cell_map = reader.get_chisel_barcodes(os.path.join(loc_dir, 'chisel_out/barcodes.info.tsv'))

    df = pd.read_csv(os.path.join(loc_dir, 'chisel_out/calls/calls.tsv'), sep='\t')
    df = df.drop(['END', 'NORM_COUNT', 'COUNT', 'A_COUNT', 'B_COUNT', 'CN_STATE'], axis=1)
    df = df.rename(columns = dict(zip(df.columns, ['chrom', 'start', 'cell', 'RDR', 'BAF', 'clust'])))
    df['cell'] = df['cell'].apply(lambda x: cell_map[x])
    df['BAF'] = df['BAF'].apply(lambda x: min(x, 1-x))
    df['start'] = df['start'] + 1

    grouped = df.groupby('clust')
    means = grouped[['RDR', 'BAF']].mean().values
    covs = np.array([group[['RDR', 'BAF']].cov().values for name, group in grouped])
    cluster_sizes = grouped.size().values
    total_size = cluster_sizes.sum()
    weights = cluster_sizes / total_size

    df['id'] = df['chrom'] + '-' + df['start'].astype(str)
    RDR_df = df.pivot(index='cell', columns='id', values='RDR')
    BAF_df = df.pivot(index='cell', columns='id', values='BAF')
    clust = df.pivot(index='cell', columns='id', values='clust')

    RDR_df = RDR_df[order]
    RDR_df = RDR_df.rename(columns = {x: i for i,x in enumerate(order)})
    RDR_df = RDR_df.reindex(cell_names)
    BAF_df = BAF_df[order]
    BAF_df = BAF_df.rename(columns = {x: i for i,x in enumerate(order)})
    BAF_df = BAF_df.reindex(cell_names)
    clust = clust[order]
    clust = clust.reindex(cell_names)

    binID_to_cell = {}
    clust_flat = clust.to_numpy().flatten()
    for i, cell in enumerate(cell_names):
        for j in range(nbins):
            binID_to_cell[i*nbins + j] = [cell, j]

    return RDR_df, BAF_df, clust_flat, means, covs, weights, binID_to_cell

def get_local_bkpts(loc_dir, bin_coords, num_bin_per_chrom):
    chrom_names = list(set([x[0] for x in bin_coords]))

    local_filepath = os.path.join(loc_dir, 'local_bkpts.tsv')
    readcount_filepath = os.path.join(loc_dir, 'readcounts.tsv')
    local_bkpts = parse_local_bkpts(local_filepath, readcount_filepath, num_bin_per_chrom)
    return local_bkpts

def get_ensemble_bkpts(cell_names, clust, local_bkpts, hfilt, binID_to_cell, nbins):
    lower_bkpts = get_bkpts(cell_names, binID_to_cell, clust, nbins)

    if hfilt > 1:
        upper_bkpts = filter_breakpoint(lower_bkpts, clust, cell_names, hfilt)
    else:
        upper_bkpts = lower_bkpts

    ensemble_bkpts = ensemble_cluster_by_frequency_filter(local_bkpts, lower_bkpts, upper_bkpts)

    for cell in cell_names:
        ensemble_bkpts[cell] = ensemble_bkpts[cell].astype(int)
    return ensemble_bkpts

def merge_static_custom(clust, means, covs, weights, t):
    
    merges = []
    best_p = get_best_pair(means, t)
    while best_p:
        c1, c2 = best_p
        merges.append((means[c1], means[c2]))
        alg_weight = weights[c1] + weights[c2]
        alg_mean = (weights[c1]*means[c1] + weights[c2]*means[c2]) / alg_weight
        alg_cov = (((weights[c1]**2)*covs[c1] + (weights[c2]**2)*covs[c2] + np.outer((weights[c1]*means[c1] + weights[c2]*means[c2]), (weights[c1]*means[c1] + weights[c2]*means[c2]).T)) / (alg_weight**2)) - np.outer(alg_mean, alg_mean.T)

        means = np.vstack((means[:c2], means[c2+1:]))
        covs = np.vstack((covs[:c2], covs[c2+1:]))
        weights = np.delete(weights, c2)

        means[c1] = alg_mean
        covs[c1] = alg_cov
        weights[c1] = alg_weight

        for i in range(len(clust)):
            if clust[i] == c2:
                clust[i] = c1
            elif clust[i] > c2:
                clust[i] -= 1
        
        best_p = get_best_pair(means, t)
    
    return clust, means, covs, weights, merges

def CN_call(out_dir, cell_names, clust, ensemble_bkpts, RDR_df, BAF_df, means, covs, weights, bin_coords, num_bin_per_chrom, max_wgd, max_CN):
    new_clust = get_segs(cell_names, clust, ensemble_bkpts, RDR_df, BAF_df, means, num_bin_per_chrom, 1)
    best_CNs, chosen = get_cluster_CNs(cell_names, [], new_clust, RDR_df, BAF_df, means, covs, weights, max_wgd, max_CN, 'unweighted')

    bin_coords = [bin_coords[i] + (i,) for i in range(len(bin_coords))]
    allele_CNs = allele_caller(cell_names, best_CNs, new_clust)

    flat_df = {}
    n_bins = len(bin_coords)
    for cell in cell_names:
        cell_CNs = allele_CNs[cell]
        flat_df[cell] = [f'{cell_CNs[i][0]},{cell_CNs[i][1]}' for i in range(n_bins)]
    flat_df = pd.DataFrame.from_dict(flat_df, orient='index')
    flat_df.index.name = 'cell'
    flat_df.to_csv(os.path.join(out_dir, 'calls_flat.tsv'), sep='\t')

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in-dir', type=str, required=True, help='Path to output directory of SEACON.')
    parser.add_argument('-c', '--chis-dir', type=str, required=True, help='Path to output directory of CHISEL.')
    parser.add_argument('-o', '--out-dir', type=str, default='chisel_heuristic', help='Output directory.')
    parser.add_argument('--upper-filter', type=int, default=5, help='Upper bound for breakpoint filter.')
    parser.add_argument('--merge-tol', type=float, default=0.1, help='Merge tolerance.')
    parser.add_argument('--max-wgd', type=int, default=1, help='Maximum number of WGDs to consider.')
    parser.add_argument('--max-CN', type=int, default=10, help='Maximum CN to consider.')
    arguments = parser.parse_args()
    return arguments

def main():
    args = parse_args()
    in_dir = args.in_dir
    out_dir = args.out_dir
    hfilt = args.upper_filter
    merge_tol = args.merge_tol
    max_wgd = args.max_wgd
    max_CN = args.max_CN

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    regions = reader.read_region_file(os.path.join(in_dir, 'filtered_regions.bed'))
    reader.write_region_file(out_dir, 'filtered_regions.bed', regions)

    cell_names = reader.read_cell_names(os.path.join(in_dir, 'cells.txt'))
    bin_coords = index_bins(os.path.join(in_dir, 'filtered_regions.bed'))
    num_bin_per_chrom = get_bins_per_chrom(bin_coords)
    nbins = len(bin_coords)
    RDR_df, BAF_df, clust, means, covs, weights, binID_to_cell = get_info(in_dir, bin_coords, cell_names)
    clust, means, covs, weights, merges = merge_static_custom(clust, means, covs, weights, merge_tol)

    local_bkpts = get_local_bkpts(in_dir, bin_coords, num_bin_per_chrom)
    ensemble_bkpts = get_ensemble_bkpts(cell_names, clust, local_bkpts, hfilt, binID_to_cell, nbins)

    CN_call(out_dir, cell_names, clust, ensemble_bkpts, RDR_df, BAF_df, means, covs, weights, bin_coords, num_bin_per_chrom, max_wgd, max_CN)

if __name__ == '__main__':
    main()
