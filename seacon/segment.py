import os
import pickle
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

from seacon.local import local_segmentation_generator
    
def euc(a, b):
    t = a - b
    d = np.sqrt(np.dot(t.T, t))
    return d

def get_bkpts(cell_names, binID_to_cell, cluster_assignments, nbins):
    """
    Args:
    - cell_names: List of cell names.
    - binID_to_cell: Dictionary the cell and cell-specific location of flattened bin ids.
    - cluster_assignments: Fully flattened list of cluster assignments.

    Returns:
    - bkpts: A dictionary of cells with cluster assignments.
    """
    bkpts = {cell: [] for cell in cell_names}
    cur_cell, cur_clust = binID_to_cell[0][0], cluster_assignments[0]
    for i,c in enumerate(cluster_assignments[1:]):
        if binID_to_cell[i][0] != cur_cell:
            cur_cell = binID_to_cell[i][0]
            cur_clust = c
        else:
            if c != cur_clust:
                if binID_to_cell[i][1]+1 != nbins:
                    bkpts[cur_cell].append(binID_to_cell[i][1]+1)
                cur_clust = c
    return bkpts

# Both pred_bkpts and pred_clust are dictionaries where keys are cell names and values correspond to flattened bins
def filter_breakpoint(pred_bkpts, pred_clust, cell_names, filt):
    num_bins = int(len(pred_clust)/len(cell_names))
    pred_bkpt_counts = {b: {} for b in range(num_bins)}
    for b in range(num_bins):
        for i,cell in enumerate(cell_names):
            if b in pred_bkpts[cell]:
                # Ensures breakpoints are counted as similar if they have the same start/end clusters
                key = str(pred_clust[(i*num_bins) + b - 1]) + '-' + str(pred_clust[(i*num_bins) + b])
                if key not in pred_bkpt_counts[b]:
                    pred_bkpt_counts[b][key] = []
                pred_bkpt_counts[b][key].append(cell)

    new_bkpts = {cell: [] for cell in cell_names}
    for b in pred_bkpt_counts:
        for k in pred_bkpt_counts[b]:
            if len(pred_bkpt_counts[b][k]) >= filt:
                for cell in pred_bkpt_counts[b][k]:
                    new_bkpts[cell].append(b)
    return new_bkpts

def get_distances(means):
    n = len(means)
    distances = [[-1 for i in range(n)] for i in range(n)]
    best_p, best_d = None, None

    for i in range(n):
        for j in range(i+1, n):
            euc_d = euc(means[i], means[j])
            distances[i][j] = distances[j][i] = euc_d #np.abs(means[i][0] - means[j][0])

            if best_d is None or euc_d < best_d:
                best_d = euc_d
                best_p = (i,j)

    return best_p, distances


def get_best_pair(means, t):
    n = len(means)
    best_p = None
    best_d = None

    for i in range(n):
        for j in range(i+1, n):
            euc_d = euc(means[i], means[j])
            if euc_d <= t:
                if best_d is None or euc_d < best_d:
                    best_p = (i,j)
                    best_d = euc_d

    return best_p


def get_gmm_model(f, min_component, max_component, init_params='kmeans', max_iter=500, n_init=1):
    models = {}
    for i in range(min_component, max_component+1):
        print(f'components={i}')
        gmm = GaussianMixture(n_components=i, init_params=init_params, max_iter=max_iter, n_init=n_init, covariance_type='full').fit(f)
        models[i] = gmm
    
    bics = [models[i].bic(f) for i in models.keys()]
    best_k, best_bic = min([(i,x) for i,x in enumerate(bics, min_component)], key= lambda x: x[1])

    return models, best_k, best_bic

def merge_relative(f, gmm, t):
    means, covs, weights = gmm.means_, gmm.covariances_, gmm.weights_
    clust = gmm.predict(f)
    p, distances = get_distances(means)
    flatds = np.triu(distances, 1).flatten()
    thresh = t*np.mean(flatds[flatds != -1])

    pairs = []
    if distances[p[0]][p[1]] <= thresh:
        pairs.append(p)

    num_merges = 0
    while pairs:
        num_merges += 1
        c1, c2 = pairs.pop()
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

        p, distances = get_distances(means)
        if distances[p[0]][p[1]] <= thresh:
            pairs.append(p)
    
    print(f'Num merges: {num_merges}')
    return clust, means, covs, weights

def merge_static(f, gmm, t):
    means, covs, weights = gmm.means_, gmm.covariances_, gmm.weights_
    clust = gmm.predict(f)
    
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
    
    print(f'Num merges: {len(merges)}')
    return clust, means, covs, weights, merges

def global_seg(flat_data, min_component, max_component, tol, gmm_params):
    f = np.array(flat_data)

    models, best_K, best_bic = get_gmm_model(f, min_component, max_component, **gmm_params)
    best_model = models[best_K]
    print('best', best_K)
    means, covs, weights = best_model.means_, best_model.covariances_, best_model.weights_
    clust, means, covs, weights, merges = merge_static(f, best_model, t=tol)

    return clust, means, covs, weights, best_K, merges


#def ensemble_cluster_by_frequency_filter(seg_local, seg_global, seg_global_filtered_dict, upper_threshold, lower_threshold=1):
def ensemble_cluster_by_frequency_filter(seg_local, lower_seg_global, upper_seg_global):
    """
    Performs ensemble clustering by frequency filtering.

    Args:
    - seg_local: Dictionary of local segments.
    - lower_seg_global: Dictionary of filtered global segments using lower bound.
    - upper_seg_global: Dictionary of filtered global segments using upper bound.

    Returns:
    - ensembled_segmentation: Dictionary of ensembled breakpoints.
    """
    ensembled_segmentation = {}
    for index, (cell, _) in enumerate(seg_local.items()):
        current_seg_local = np.array(seg_local[cell])
        current_upper_seg_global = np.array(upper_seg_global[cell])
        current_lower_seg_global = np.array(lower_seg_global[cell])

        breakpoints_in_question = np.sort(np.setdiff1d(current_seg_local, current_upper_seg_global))
        final_breakpoints = current_upper_seg_global

        for breakpoint in breakpoints_in_question:
            if breakpoint in current_lower_seg_global:
                final_breakpoints = np.append(final_breakpoints, breakpoint)
        final_breakpoints = np.sort(final_breakpoints)
        ensembled_segmentation[cell] = final_breakpoints

    return ensembled_segmentation

def parse_local_bkpts(local_filepath, readcount_filepath, num_bin_per_chrom):
    cell_names = []
    with open(readcount_filepath) as f:
        next(f)
        for line in f:
            cell_names.append(line.strip().split('\t')[0])
    
    bkpts = []
    with open(local_filepath) as f:
        for line in f:
            info = line.strip().split('\t')
            cur_bkpt = [int(x) for x in info[1:-1]]
            bkpts.append(cur_bkpt)
    local_bkpts = dict(zip(cell_names, bkpts))
    return local_bkpts


def segmentation(args, gmm_params={}, cbs_params={}):
    out_dir = args['out_dir']
    flat_data = args['flat_data']
    binID_to_cell = args['binID_to_cell']
    cell_to_binID = args['cell_to_binID']
    cell_names = args['cell_names']
    chrom_names = args['chrom_names']
    num_bin_per_chrom = args['num_bin_per_chrom']
    min_component = args['min_component']
    max_component = args['max_component']
    tol = args['tol']
    bkpt_filt_high = args['bkpt_filt_high']
    seg_type = args['seg_type']

    local_bkpts, lower_bkpts, upper_bkpts, ensemble_bkpts = None, None, None, None
    nbins = int(len(flat_data)/len(cell_names))

    if seg_type == 'local' or seg_type == 'ensemble':
        local_filepath = os.path.join(out_dir, 'local_bkpts.tsv')
        readcount_filepath = os.path.join(out_dir, 'readcounts.tsv')
        local_segmentation_generator(local_filepath, readcount_filepath, len(chrom_names), num_bin_per_chrom=num_bin_per_chrom, **cbs_params)
        local_bkpts = parse_local_bkpts(local_filepath, readcount_filepath, num_bin_per_chrom)

        if seg_type == 'local':
            return {'clust': None, 'means': None, 'covs': None, 'weights': None, 'best_K': None, 'bkpts': local_bkpts, 'assignments': None}
    
    clust, means, covs, weights, best_K, merges = global_seg(flat_data, min_component, max_component, tol, gmm_params)

    cell_clust = dict(zip(cell_names, np.split(clust, len(cell_names))))

    lower_bkpts = get_bkpts(cell_names, binID_to_cell, clust, nbins)

    if bkpt_filt_high > 1:
        upper_bkpts = filter_breakpoint(lower_bkpts, clust, cell_names, bkpt_filt_high)
    else:
        upper_bkpts = lower_bkpts

    if seg_type == 'global':
        ensemble_bkpts = upper_bkpts
    elif seg_type == 'ensemble':
        ensemble_bkpts = ensemble_cluster_by_frequency_filter(local_bkpts, lower_bkpts, upper_bkpts)

    for cell in cell_names:
        ensemble_bkpts[cell] = ensemble_bkpts[cell].astype(int)

    return {'clust': clust, 'means': means, 'covs': covs, 'weights': weights, 'best_K': best_K, 'bkpts': ensemble_bkpts, 'merges': merges}