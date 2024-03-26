import numpy as np
import pandas as pd
import os
import argparse
import datetime
import pickle
import pathlib

import reader
import normalizer as norm
from helpers import *
from segment import *
from caller import *

def parse_args():
    parser = argparse.ArgumentParser()

    # Basic stuff
    parser.add_argument('mode', type=str, help='Options are \'prep\', \'norm\', \'pseudonormal\', \'call\', and \'full\'.')
    parser.add_argument('-o', '--out-dir', type=str, default=None, help='Output directory. Default: current directory.')
    parser.add_argument('-i', '--bam-path', type=str, default=None, help='Path to directory containing bam files.')
    parser.add_argument('-r', '--reference', type=str, default=None, help='Path to reference genome.')
    parser.add_argument('-b', '--bin-size', type=int, default=None, help='Bin size in bp.')

    # Inputs and Executables
    parser.add_argument('--bedtools', type=str, default='bedtools', help='Path to bedtools executable. Default: in user $PATH.')
    parser.add_argument('--bigwig', type=str, default='bigWigAverageOverBed', help='Path to bigWigAverageOverBed executable. Default: in user $PATH.')
    parser.add_argument('--map-file', type=str, default='assets/hg38_mappability.bigWig', help='Path to mappability file in bigwig format. By default, set to assets/hg38_mappability.bigWig.')
    parser.add_argument('--baf-path', type=str, default=None, help='Path to input BAFs. If using BAFs derived by CHISEL, set path to the output directory of CHISEL. Otherwise, set path to the required tsv file.')

    # Normalization
    parser.add_argument('--no-normal', action='store_true', help='Indicates no normal cells in sample. Applies standard gc and mappability corrections instead.')
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

    parser.add_argument('-m4', action='store_true', help='Testing.')

    arguments = parser.parse_args()
    checked_args = check_args(arguments)
    return checked_args

def check_args(arguments):
    args = vars(arguments)
    if not os.path.isdir(args['out_dir']):
        os.mkdir(args['out_dir'])

    #if args['map_file'] == None:
    #    source_dir = pathlib.Path(__file__).parent.resolve()
    #    args['map_file'] = os.path.join(source_dir, 'assets/hg38_mappability.bigWig')
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

def main(args):
    save_args(args)
    if args['mode'] == 'prep' or args['mode'] == 'full':
        chrom_lens = reader.get_chrom_lens(args['reference'])
        chrom_names = list(chrom_lens.keys())
        print(chrom_names)

        regions = reader.get_bins_from_chromlens(chrom_lens, args['bin_size'])
        full_region_path = reader.write_region_file(args['out_dir'], 'full_regions.bed', regions)

        gc = reader.get_gc(args['reference'], full_region_path, args['bedtools'])
        mapp = reader.get_mapp(args['out_dir'], full_region_path, args['bigwig'], args['map_file'])

        filtered_regions = reader.filter_bins_gc_mapp(regions, gc, mapp)
        filtered_region_path = reader.write_region_file(args['out_dir'], 'filtered_regions.bed', filtered_regions)

        if os.path.isdir(args['bam_path']):
            raw_readcounts_df, cell_names = reader.get_readcounts_bamdir(args['bam_path'], filtered_region_path)
        else:
            print('Cannot find bam dir.')
            return
        print('Number of cells:', len(cell_names))
        raw_readcounts_df.to_csv(os.path.join(args['out_dir'], 'raw_readcounts.tsv'), sep='\t')  
        reader.write_cell_names(os.path.join(args['out_dir'], 'cells.txt'), cell_names)

    if args['mode'] == 'norm' or args['mode'] == 'full':
        cell_names = reader.read_cell_names(os.path.join(args['out_dir'], 'cells.txt'))
        bin_coords = index_bins(os.path.join(args['out_dir'], 'filtered_regions.bed'))
        raw_readcounts_df = reader.read_data_flat(os.path.join(args['out_dir'], 'raw_readcounts.tsv'))
        if args['no_normal']:
            print('Not implemented.')
        else:
            ginis = norm.compute_gini_samples(raw_readcounts_df)
            with open(os.path.join(args['out_dir'], 'gini.pkl'), 'wb') as f:
                pickle.dump(ginis, f)
            if args['gini_thresh']:
                normal_cells = norm.get_normal_from_gini_thresh(ginis, args['gini_thresh'])
            else:
                normal_cells = norm.get_normal_from_gini_auto(ginis, min_count = args['min_norm'])
            with open(os.path.join(args['out_dir'], 'diploid.txt'), 'w+') as f:
                for cell in normal_cells:
                    f.write(cell + '\n')

            readcounts_df = norm.normcell_correct(raw_readcounts_df, normal_cells)
            readcounts_df[readcounts_df.columns] = readcounts_df[readcounts_df.columns].astype(int)
            readcounts_df.to_csv(os.path.join(args['out_dir'], 'readcounts.tsv'), sep='\t')

    # Must have already run method in 'norm' mode, or create file in out_dir named 'diploid.txt' which has the cell names of each diploid cell
    if args['mode'] == 'pseudonormal':
        if os.path.isdir(args['bam_path']):
            norm.make_pseudonormal(args['bam_path'], args['out_dir'])
        else:
            print('Cannot find bam dir.')
            return

    if args['mode'] == 'call' or args['mode'] == 'full':
        cell_names = reader.read_cell_names(os.path.join(args['out_dir'], 'cells.txt'))
        bin_coords = index_bins(os.path.join(args['out_dir'], 'filtered_regions.bed'))
        readcounts_df = reader.read_data_flat(os.path.join(args['out_dir'], 'readcounts.tsv'))
        if os.path.exists(os.path.join(args['out_dir'], 'diploid.txt')):
            normal_cells = reader.read_cell_names(os.path.join(args['out_dir'], 'diploid.txt'))
        else:
            normal_cells = []

        if os.path.exists(os.path.join(args['out_dir'], 'RDR.tsv')) and os.path.exists(os.path.join(args['out_dir'], 'BAF.tsv')) and args['baf_path'] == None:
            RDR_df = reader.read_data_flat(os.path.join(args['out_dir'], 'RDR.tsv'))
            BAF_df = reader.read_data_flat(os.path.join(args['out_dir'], 'BAF.tsv'))
        else:
            if args['baf_path'] == None:
                print('Need to specify path to BAFs.')
                return
            else:
                if os.path.isdir(args['baf_path']):
                    cellmap = reader.get_chisel_barcodes(os.path.join(args['baf_path'], 'barcodes.info.tsv'))
                    BAF_df, alt_bin_coords = reader.collect_chisel_BAF(os.path.join(args['baf_path'], 'combo/combo.tsv'), cell_names, cellmap = cellmap)
                else:
                    BAF_df, alt_bin_coords = reader.collect_other_BAF(args['baf_path'], cell_names)

                readcounts_df, BAF_df, shared_bins, shared_cells = correct_RC_BAF(readcounts_df, BAF_df, bin_coords, alt_bin_coords, normal_cells=normal_cells)
                cell_avgs = readcounts_df.mean(axis=1)
                RDR_df = readcounts_df.div(cell_avgs, axis=0)

                readcounts_df[readcounts_df.columns] = readcounts_df[readcounts_df.columns].astype(int)
                readcounts_df.to_csv(os.path.join(args['out_dir'], 'readcounts.tsv'), sep='\t')

                filtered_region_path = reader.write_region_file(args['out_dir'], 'filtered_regions.bed', shared_bins)
                reader.write_cell_names(os.path.join(args['out_dir'], 'cells.txt'), shared_cells)
                cell_names = shared_cells
                RDR_df = RDR_df.round(5)
                BAF_df = BAF_df.round(5)
                RDR_df.to_csv(os.path.join(args['out_dir'], 'RDR.tsv'), sep='\t')
                BAF_df.to_csv(os.path.join(args['out_dir'], 'BAF.tsv'), sep='\t')

        chrom_names = list(set([x[0] for x in bin_coords]))
        chrom_names.sort(key = lambda x: int(x[3:]))
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
        with open(os.path.join(args['out_dir'], 'seg_obj.pkl'), 'wb') as f:
            pickle.dump(seg_obj, f)
    
    if args['mode'] == 'call' or args['mode'] == 'full' or args['mode'] == 'infer':

        if args['mode'] == 'infer':
            RDR_df = reader.read_data_flat(os.path.join(args['out_dir'], 'RDR.tsv'))
            BAF_df = reader.read_data_flat(os.path.join(args['out_dir'], 'BAF.tsv'))

            cell_names = list(RDR_df.index)
            bin_coords = index_bins(os.path.join(args['out_dir'], 'filtered_regions.bed'))
            num_bin_per_chrom = get_bins_per_chrom(bin_coords)
            chrom_names = list(set([x[0] for x in bin_coords]))
            chrom_names.sort(key = lambda x: int(x[3:]))
            flat_data, binID_to_cell, cell_bin_to_binID = combine_flat_data(RDR_df, BAF_df)

            with open(os.path.join(args['out_dir'], 'seg_obj.pkl'), 'rb') as f:
                seg_obj = pickle.load(f)

        ensemble_bkpts = seg_obj['bkpts']
        cell_names = list(ensemble_bkpts.keys())
        clust = seg_obj['clust']
        means = seg_obj['means']
        covs = seg_obj['covs']
        weights = seg_obj['weights']

        if os.path.exists(os.path.join(args['out_dir'], 'diploid.txt')):
            normal_cells = reader.read_cell_names(os.path.join(args['out_dir'], 'diploid.txt'))
        else:
            normal_cells = []



        new_clust = get_segs(cell_names, clust, ensemble_bkpts, RDR_df, BAF_df, means, num_bin_per_chrom, args['min_seg_length'])

        best_CNs, chosen = get_cluster_CNs(cell_names, normal_cells, new_clust, RDR_df, BAF_df, means, covs, args['max_wgd'], args['max_CN'], args['ploidy_RDR'])

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
            frames.append(pd.DataFrame(cell_frames, columns=['cell', 'chrom', 'start', 'end', 'hapCN']))
        df = pd.concat(frames)
        df.to_csv(os.path.join(args['out_dir'], f'calls.tsv'), sep='\t', index=False)



if __name__ == '__main__':
    args = parse_args()
    main(args)