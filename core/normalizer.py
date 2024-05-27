import numpy as np
import pandas as pd
import pysam
from scipy.signal import find_peaks

import os

# Assume we have already removed bad bins and reformated data to be a list CNP for each cell. Normalizes read counts based on mappability.
def map_correct(readcount_df, mapp):
    good_bins = list(readcount_df.columns)
    for b in good_bins:
        readcount_df[b] = readcount_df[b] / mapp[b]
    return readcount_df

# Normalizes read counts based on GC content bias.
def GC_correct(readcount_df, gc):
    # Arrange bins by GC content
    GC_percentiles = np.percentile(gc, range(0, 105, 5))
    GC_bins = [[] for i in range(len(GC_percentiles)-1)]
    for i in range(len(gc)):
        for j in range(len(GC_percentiles)-1):
            if gc[i] >= GC_percentiles[j] and gc[i] < GC_percentiles[j+1]:
                GC_bins[j].append(i)
                break
            elif gc[i] == GC_percentiles[-1]:
                GC_bins[-1].append(i)

    # Correct for each cell independently
    temp_dfs = []
    for group in GC_bins:
        mean_vals = readcount_df[group].mean(axis=1)
        temp_df = pd.DataFrame({j: mean_vals for j in group}, index=readcount_df.index)
        temp_dfs.append(temp_df)
        
    GC_bin_avg = pd.concat(temp_dfs, axis=1)
    GC_bin_avg = GC_bin_avg.sort_index(axis=1)
    cell_means = readcount_df.mean(axis=1)
    readcount_df = readcount_df.mul(cell_means, axis=0)
    readcount_df = readcount_df.div(GC_bin_avg)
            
    return readcount_df

# Normalizes based on reads accumulated from normal cells
def normcell_correct(readcount_df, normal_cells):
    norm_rc = readcount_df.filter(items=normal_cells, axis=0)
    norm_means = norm_rc.apply(np.mean, axis=1)
    # First divide by bin-lambda
    for b in range(len(readcount_df.columns)):
        norm_bin_vals = norm_rc[b] / norm_means
        bin_lambda = np.mean(norm_bin_vals)
        readcount_df[b] = readcount_df[b] / bin_lambda
    return readcount_df

def bulknorm_correct(readcount_df, norm_counts):
    norm_mean = np.mean(norm_counts)
    bin_lambdas = norm_counts / norm_mean
    for b in range(len(readcount_df.columns)):
        readcount_df[b] = readcount_df[b] / y[b]
    return readcount_df

# Given array of read counts, computes gini coefficient
def compute_gini(cell_counts):
    total = 0
    for i, x in enumerate(cell_counts[:-1], 1):
        total += np.sum(np.abs(x - cell_counts[i:]))
    gini = total / (len(cell_counts)**2 * np.mean(cell_counts))
    return gini

# Creates series of gini coefficients for each cell
def compute_gini_samples(readcounts_df):
    gini = readcounts_df.apply(compute_gini, axis=1)
    return gini

# Returns list of normal cell names if gini below threshold
def get_normal_from_gini_thresh(gini, thresh):
    subset = gini[gini < thresh]
    normal_cell = list(subset.index)
    return normal_cell

# Automatically selects normal cells as first peak of gini coefficients
def get_normal_from_gini_auto(gini, num_boxes = None, min_count = 10, max_gini = 0.3):
    gini_dict = gini.to_dict()
    gini_flattened = list(gini_dict.items())
    gini_flattened.sort(key = lambda x: x[1])
    [cell_names, g_vals] = map(list, zip(*gini_flattened))

    minVal, maxVal = min(g_vals), max(g_vals)
    num_cells = len(cell_names)

    if num_boxes == None:
        num_boxes = 100
    interval_len = 1 / num_boxes
    intervals = [interval_len*i for i in range(1, num_boxes+1)]
    boxes = [[] for i in range(num_boxes)]

    cur_idx, cur_box = 0, 0
    while cur_idx < num_cells:
        while cur_box < len(intervals) and g_vals[cur_idx] > intervals[cur_box]:
            cur_box += 1
        if cur_box == len(intervals):
            boxes[cur_box].extend([i for i in range(cur_idx, num_cells)])
            break
        else:
            boxes[cur_box].append(cur_idx)
        cur_idx += 1
    
    box_counts = [len(x) for x in boxes]
    min_height = max(round(num_cells*0.05), 1)
    peaks = find_peaks(box_counts, height=min_height)
    if len(peaks[0]) < 2:
         peaks = find_peaks(box_counts)
    if len(peaks[0]) == 0:
        print('Error')
        return

    if len(peaks[0]) == 1:
        end_idx = peaks[0][0] + 1
        num_normal = sum(box_counts[:end_idx])
        while num_normal < min_count:
            end_idx += 1
            num_normal = sum(box_counts[:end_idx])
    else:
        i = 1
        end_idx = int(np.floor((peaks[0][i-1] + peaks[0][i])/2))
        num_normal = sum(box_counts[:end_idx])
        while i+1 < len(peaks[0]) and num_normal <= min_count: 
            i += 1
            end_idx = np.floor((peaks[0][i-1] + peaks[0][i])/2)
            num_normal = sum(box_counts[:end_idx])

    normal_cells = []
    for i in range(end_idx):
        for j in boxes[i]:
            if g_vals[j] <= max_gini:
                normal_cells.append(cell_names[j])

    return normal_cells

def make_pseudonormal(bam_dir, out_dir, normal_cells):
    temp_paths = os.path.join(out_dir, 'norm_paths.txt')
    with open(temp_paths, 'w+') as f:
        for cell in normal_cells:
            cur_path = os.path.join(bam_dir, cell + '.bam')
            f.write(f'{cur_path}\n')
    
    pseudo_path = os.path.join(out_dir, 'pseudonormal.bam')
    pysam.merge('-o', pseudo_path, '-b', temp_paths)
    pysam.index(pseudo_path)
    os.remove(temp_paths)