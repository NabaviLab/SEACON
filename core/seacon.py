#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os
import pathlib
import argparse
import datetime
import pickle

import core.reader as reader
import core.normalizer as norm
from core.helpers import *
from core.segment import *
from core.caller import *
from core.baf import chisel_baf_helper, convert_baf, call_chisel_combo, call_chisel_caller

def parse_args():
    parser = argparse.ArgumentParser()

    # Basics, Inputs and Executables
    parser.add_argument('mode', type=str, help='Options are \'prep_readcount\', \'pseudonormal\', \'prep_baf\', \'call\', and \'full\'.')
    parser.add_argument('-o', '--out-dir', type=str, default=None, help='Output directory. Default: current directory.')
    parser.add_argument('-i', '--bam-path', type=str, default=None, help='Path to directory containing bam files.')
    parser.add_argument('-r', '--reference', type=str, default=None, help='Path to reference genome.')
    parser.add_argument('--chrom_names', type=str, default=None, help='Path to file containing names of chromosomes, one per line.')
    parser.add_argument('-b', '--bin-size', type=int, default=None, help='Bin size in bp.')
    parser.add_argument('--bedtools', type=str, default='bedtools', help='Path to bedtools executable. Default: in user $PATH.')
    parser.add_argument('--bigwig', type=str, default='bigWigAverageOverBed', help='Path to bigWigAverageOverBed executable. Default: in user $PATH.')
    parser.add_argument('--map-file', type=str, default='assets/hg38_mappability.bigWig', help='Path to mappability file in bigwig format. By default, set to assets/hg38_mappability.bigWig.')
    parser.add_argument('--precomputed-baf', type=str, default=None, help='Path to file containing precomputed BAFs. See instructions for file format.')
    parser.add_argument('--vcf', type=str, default=None, help='Path to input vcf file containing phased SNPs.')
    parser.add_argument('--block-size', type=int, default=100, help='Size of haplotype blocks in kb (default 100).')
    parser.add_argument('-P', '--num-processors', type=int, default=1, help='Number of parallel processors to use. Defaults to 1.')
    
    # Normalization
    parser.add_argument('--no-normal', action='store_true', help='Indicates no normal cells in sample. Applies standard gc and mappability corrections instead.')
    parser.add_argument('--normal-bam', type=str, default=None, help='Path to matched-normal bam file.')
    parser.add_argument('--gini-thresh', type=float, default=None, help='Strict gini threshold to use for normal cell identification. Default: determined automatically.')
    parser.add_argument('--min-norm', type=int, default=10, help='Minimum number of normal cells to use.')

    # Segmentation and Calling
    parser.add_argument('--min-component', type=int, default=8, help='Minimum number of GMM components to consider.')
    parser.add_argument('--max-component', type=int, default=25, help='Maximum number of GMM components to consider.')
    parser.add_argument('--max-iter', type=int, default=1000, help='Max number of GMM iterations.')
    parser.add_argument('--tolerance', type=float, default=0.2, help='Tolerance for combining GMM components.')
    parser.add_argument('--upper-filter', type=int, default=5, help='Upper bound for breakpoint filter.')
    parser.add_argument('--max-wgd', type=int, default=2, help='Maximum number of WGDs to consider.')
    parser.add_argument('--max-CN', type=int, default=10, help='Maximum CN to consider.')
    parser.add_argument('--ploidy-RDR', action='store_true', help='Select ploidy from RDR instead of combined RDR/BAF. Use if BAF estimates are highly noisy.')
    parser.add_argument('--num-permutation', type=int, default=1000, help='Number of permutations in CBS algorithm.')
    parser.add_argument('--alpha', type=float, default=0.05, help='Alpha value for CBS algorithm.')
    parser.add_argument('--min-width', type=int, default=3, help='Minimum width for CBS algorithm.')
    parser.add_argument('--k-max', type=int, default=10, help='Maximum number of change-points for CBS algorithm.')
    parser.add_argument('--min-seg-length', type=int, default=1, help='Minimum number of bins in a continuous segment.')

    arguments = parser.parse_args()
    checked_args = check_args(arguments)
    return checked_args

def check_args(arguments):
    args = vars(arguments)
    if not os.path.isdir(args['out_dir']):
        os.mkdir(args['out_dir'])
    return args

def save_args(args):
    if os.path.isfile(os.path.join(args['out_dir'], 'log.txt')):
        with open(os.path.join(args['out_dir'], 'log.txt'), 'a') as f:
            cur_time = datetime.datetime.now()
            f.write(f'{cur_time}\n')
            for i,v in args.items():
                f.write(f'{i}:\t{v}\n')
            f.write('\n')
    else:
        with open(os.path.join(args['out_dir'], 'log.txt'), 'w+') as f:
            cur_time = datetime.datetime.now()
            f.write(f'{cur_time}\n')
            for i,v in args.items():
                f.write(f'{i}:\t{v}\n')
            f.write('\n')

def prep_readcount(args):
    chrom_names = reader.get_chrom_names(args['chrom_names'])
    chrom_lens = reader.get_chrom_lens(args['reference'], chrom_names)
    regions = reader.get_bins_from_chromlens(chrom_lens, args['bin_size'])
    full_region_path = reader.write_region_file(args['out_dir'], 'full_regions.bed', regions)

    gc = reader.get_gc(args['reference'], full_region_path, args['bedtools'])
    mapp = reader.get_mapp(args['out_dir'], full_region_path, args['bigwig'], args['map_file'])

    filtered_regions, filtered_stats = reader.filter_bins_gc_mapp(regions, gc, mapp)
    reader.write_region_file(args['out_dir'], 'filtered_regions.bed', filtered_regions, stats=filtered_stats)
    bin_coords = index_bins(os.path.join(args['out_dir'], 'filtered_regions.bed'))

    raw_readcounts_df, cell_names = reader.get_readcounts_bamdir(args['bam_path'], filtered_regions, num_processors=args['num_processors'])
    print('Number of cells:', len(cell_names))

    raw_readcounts_df.to_csv(os.path.join(args['out_dir'], 'raw_readcounts.tsv'), sep='\t')  
    reader.write_cell_names(os.path.join(args['out_dir'], 'cells.txt'), cell_names)

    if args['no_normal']:
        gc, mapp = reader.read_stats(filtered_regions)
        readcounts_df = norm.map_correct(raw_readcounts_df, mapp)
        readcounts_df = norm.GC_correct(readcounts_df, gc)
    else:
        if args['normal_bam']:
            norm_counts = reader.get_readcounts_bamfile(args['normal_bam'], filtered_regions)
            readcounts_df = norm.bulknorm_correct(raw_readcounts_df, norm_counts)
        else:
            ginis = norm.compute_gini_samples(raw_readcounts_df)
            #with open(os.path.join(args['out_dir'], 'gini.pkl'), 'wb') as f:
            #    pickle.dump(ginis, f)
            if args['gini_thresh']:
                normal_cells = norm.get_normal_from_gini_thresh(ginis, args['gini_thresh'])
            else:
                normal_cells = norm.get_normal_from_gini_auto(ginis, min_count = args['min_norm'])
            reader.write_cell_names(os.path.join(args['out_dir'], 'diploid.txt'), normal_cells)
            readcounts_df = norm.normcell_correct(raw_readcounts_df, normal_cells)

    cell_avgs = readcounts_df.mean(axis=1)
    RDR_df = readcounts_df.div(cell_avgs, axis=0)
    RDR_df = RDR_df.round(5)
    RDR_df.to_csv(os.path.join(args['out_dir'], 'RDR.tsv'), sep='\t')
    readcounts_df[readcounts_df.columns] = readcounts_df[readcounts_df.columns].astype(int)
    readcounts_df.to_csv(os.path.join(args['out_dir'], 'readcounts.tsv'), sep='\t')

def pseudonormal(args):
    normal_cells = reader.read_cell_names(os.path.join(args['out_dir'], 'diploid.txt'))
    norm.make_pseudonormal(args['bam_path'], args['out_dir'], normal_cells)

def prep_baf(args):
    cell_names = reader.read_cell_names(os.path.join(args['out_dir'], 'cells.txt'))
    chrom_names = reader.get_chrom_names(args['chrom_names'])
    bin_coords = index_bins(os.path.join(args['out_dir'], 'filtered_regions.bed'))
    readcounts_df = reader.read_data_flat(os.path.join(args['out_dir'], 'readcounts.tsv'))
    if args['precomputed_baf']:
        BAF_df, alt_bin_coords = reader.collect_other_BAF(args['precomputed_baf'], cell_names)
        readcounts_df, BAF_df, shared_bins, shared_cells = correct_RC_BAF(readcounts_df, BAF_df, bin_coords, alt_bin_coords)
        BAF_df = BAF_df.round(5)
        BAF_df.to_csv(os.path.join(args['out_dir'], 'BAF.tsv'), sep='\t')
        if len(shared_bins) < len(bin_coords):
            reader.write_region_file(args['out_dir'], 'filtered_regions.bed', shared_bins, stats=filtered_stats)
            reader.write_cell_names(os.path.join(args['out_dir'], 'cells.txt'), shared_cells)
            cell_avgs = readcounts_df.mean(axis=1)
            RDR_df = readcounts_df.div(cell_avgs, axis=0)
            RDR_df = RDR_df.round(5)
            RDR_df.to_csv(os.path.join(args['out_dir'], 'RDR.tsv'), sep='\t')
    else:
        temp_chis_dir = os.path.join(args['out_dir'], 'chisel_out')
        if not os.path.isdir(temp_chis_dir):
            os.mkdir(temp_chis_dir)
        RDR_df = reader.read_data_flat(os.path.join(args['out_dir'], 'RDR.tsv'))
        if args['normal_bam']:
            norm_counts = reader.get_readcounts_bamfile(args['normal_bam'], filtered_regions)
        else:
            normal_cells = reader.read_cell_names(os.path.join(args['out_dir'], 'diploid.txt'))
            norm_cell_rcs = readcounts_df[readcounts_df.index.isin(normal_cells)]
            norm_counts = []
            for i in range(len(bin_coords)):
                norm_counts.append(np.sum(norm_cell_rcs[i]))

        if 'CHIS_SRC' in os.environ:
            chisel_src = os.environ['CHIS_SRC']
        else:
            root_dir = pathlib.Path(__file__).parent.parent.resolve()
            chisel_src = os.path.join(root_dir, 'chisel/src/chisel')
            if not os.path.exists(chisel_src):
                print('Cannot find chisel submodule.')
                return

        chisel_baf_helper(cell_names, chrom_names, bin_coords, norm_counts, readcounts_df, RDR_df, args['bam_path'], args['vcf'], args['out_dir'], temp_chis_dir, args['num_processors'])
        call_chisel_combo(chisel_src, temp_chis_dir, args['block_size'], args['num_processors'])
        call_chisel_caller(chisel_src, temp_chis_dir, args['max_CN'], args['num_processors'])
        convert_baf(args['out_dir'], temp_chis_dir, cell_names, bin_coords)

def call(args):
    cell_names = reader.read_cell_names(os.path.join(args['out_dir'], 'cells.txt'))
    bin_coords = index_bins(os.path.join(args['out_dir'], 'filtered_regions.bed'))
    readcounts_df = reader.read_data_flat(os.path.join(args['out_dir'], 'readcounts.tsv'))

    if os.path.exists(os.path.join(args['out_dir'], 'RDR.tsv')) and os.path.exists(os.path.join(args['out_dir'], 'BAF.tsv')):
        RDR_df = reader.read_data_flat(os.path.join(args['out_dir'], 'RDR.tsv'))
        BAF_df = reader.read_data_flat(os.path.join(args['out_dir'], 'BAF.tsv'))
    else:
        print('RDRs and BAFs not found. Make sure the previous steps are run.')
        return

    chrom_names = list(set([x[0] for x in bin_coords]))
    flat_data, binID_to_cell, cell_bin_to_binID = combine_flat_data(RDR_df, BAF_df)
    num_bin_per_chrom = get_bins_per_chrom(bin_coords)

    seg_params = {}
    seg_params['out_dir'] = args['out_dir']
    seg_params['flat_data'] = flat_data
    seg_params['binID_to_cell'] = binID_to_cell
    seg_params['cell_to_binID'] = cell_bin_to_binID
    seg_params['cell_names'] = cell_names
    seg_params['chrom_names'] = chrom_names
    seg_params['num_bin_per_chrom'] = num_bin_per_chrom
    seg_params['min_component'] = args['min_component']
    seg_params['max_component'] = args['max_component']
    seg_params['tol'] = args['tolerance']
    seg_params['bkpt_filt_high'] = args['upper_filter']
    seg_params['seg_type'] = 'ensemble'

    gmm_params = {
        'max_iter': args['max_iter']
    }

    cbs_params = {
        'num_permutation': args['num_permutation'],
        'alpha': args['alpha'],
        'min_width': args['min_width'], 
        'k_max': args['k_max']
    }

    seg_obj = segmentation(seg_params, gmm_params=gmm_params, cbs_params=cbs_params)
    print('Best k: ', seg_obj['best_K'], len(seg_obj['means']))
    #with open(os.path.join(args['out_dir'], 'seg_obj.pkl'), 'wb') as f:
    #    pickle.dump(seg_obj, f)

    #with open(os.path.join(args['out_dir'], 'seg_obj.pkl'), 'rb') as f:
    #    seg_obj = pickle.load(f)

    ensemble_bkpts = seg_obj['bkpts']
    cell_names = list(ensemble_bkpts.keys())
    clust = seg_obj['clust']
    means = seg_obj['means']
    covs = seg_obj['covs']
    weights = seg_obj['weights']

    new_clust = get_segs(cell_names, clust, ensemble_bkpts, RDR_df, BAF_df, means, num_bin_per_chrom, args['min_seg_length'])
    best_CNs, chosen = get_cluster_CNs(cell_names, new_clust, RDR_df, BAF_df, means, covs, args['max_wgd'], args['max_CN'], args['ploidy_RDR'])

    bin_coords = [bin_coords[i] + (i,) for i in range(len(bin_coords))]
    allele_CNs = allele_caller(cell_names, best_CNs, new_clust)

    flat_df = {}
    n_bins = len(bin_coords)
    for cell in cell_names:
        cell_CNs = allele_CNs[cell]
        flat_df[cell] = [f'{cell_CNs[i][0]},{cell_CNs[i][1]}' for i in range(n_bins)]
    flat_df = pd.DataFrame.from_dict(flat_df, orient='index')
    flat_df.index.name = 'cell'
    flat_df.to_csv(os.path.join(args['out_dir'], f'calls_flat.tsv'), sep='\t')

    frames = []
    for cell in cell_names:
        cell_CNs = allele_CNs[cell]
        cell_frames = [[cell, x[0], x[1], x[2], f'{cell_CNs[x[3]][0]},{cell_CNs[x[3]][1]}'] for x in bin_coords]
        frames.append(pd.DataFrame(cell_frames, columns=['cell', 'chrom', 'start', 'end', 'CN']))
    df = pd.concat(frames)
    df.to_csv(os.path.join(args['out_dir'], f'calls.tsv'), sep='\t', index=False)

def main():
    args = parse_args()
    save_args(args)
    if args['mode'] == 'prep_readcount' or args['mode'] == 'full':
        if check_prep_readcount(args) == 0:
            prep_readcount(args)

    if args['mode'] == 'pseudonormal':
        pseudonormal(args)

    if args['mode'] == 'prep_baf' or args['mode'] == 'full':
        if check_prep_baf(args) == 0:
            prep_baf(args)

    if args['mode'] == 'call' or args['mode'] == 'full':
        call(args)

if __name__ == '__main__':
    main()
