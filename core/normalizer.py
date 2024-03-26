import numpy as np
import pandas as pd
import pysam
from scipy.signal import find_peaks

import os

# Assume we have already removed bad bins and reformated data to be a list CNP for each cell. Normalizes read counts based on mappability.
def map_correct(readcount_df, map_scores):
    good_bins = list(readcount_df.columns)
    for b in good_bins:
        readcount_df[b] = readcount_df[b] / map_scores[b][0]
    return readcount_df

# Normalizes read counts based on GC content bias.
def GC_correct(readcount_df, map_scores, cell_names):
    # Arrange bins by GC content
    num_intervals = 10 if len(map_scores) < 1000 else 1
    GC_percentiles = np.percentile([k[2] for k in map_scores], range(0, 100 + num_intervals, num_intervals))
    GC_bins = [[] for i in range(len(GC_percentiles)-1)]
    for i in range(len(map_scores)):
        for j in range(0, num_intervals):
            if map_scores[i][2] >= GC_percentiles[j] and map_scores[i][2] <= GC_percentiles[j+1]:
                GC_bins[j].append(i)
                break

    # Correct for each cell independently
    for cell in cell_names:
        cell_row = readcount_df.loc[cell]
        avg_map = np.average(cell_row)
        GC_bin_avg = {}
        for group in GC_bins:
            group_avg_rc = sum(cell_row[k] for k in group) / len(group)
            for k in group:
                GC_bin_avg[k] = group_avg_rc
        
        readcount_df.loc[cell] *= avg_map
        for b in range(len(readcount_df.columns)):
            readcount_df[b].loc[cell] = readcount_df[b].loc[cell] / GC_bin_avg[b]
            
    return readcount_df

# Normalizes based on reads accumulated from normal cells
def normcell_correct(readcount_df, normal_cells):
    norm_rc = readcount_df.filter(items=normal_cells, axis=0)
    # First divide by bin-lambda
    bin_lambdas = {}
    for b in range(len(readcount_df.columns)):
        norm_bin_vals = norm_rc[b] / norm_rc.apply(np.mean, axis=1)
        bin_lambda = np.sum(norm_bin_vals) / len(norm_bin_vals)
        readcount_df[b] = readcount_df[b] / bin_lambda
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

def make_pseudonormal(bam_dir, out_dir):
    normal_cells = []
    with open(os.path.join(out_dir, 'diploid.txt')) as f:
        for line in f:
            normal_cells.append(line.strip())

    temp_paths = os.path.join(out_dir, 'norm_paths.txt')
    with open(temp_paths, 'w+') as f:
        for cell in normal_cells:
            cur_path = os.path.join(bam_dir, cell + '.bam')
            f.write(f'{cur_path}\n')
    
    pseudo_path = os.path.join(out_dir, 'pseudonormal.bam')
    pysam.merge('-o', pseudo_path, '-b', temp_paths)
    pysam.index(pseudo_path)
    os.remove(temp_paths)
