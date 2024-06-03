import numpy as np
from itertools import combinations_with_replacement
from scipy.stats import multivariate_normal
import pickle

def euc(a, b):
    t = a - b
    d = np.sqrt(np.dot(t.T, t))
    return d

# cell-specific
def gen_candidates(clust, means, weights, max_wgd, RDR):
    # Filters out small components with weight < 0.1 and RDR < 0.25.
    m = [(i,m) for i,m in enumerate(means) if weights[i] >= 0.1 and means[i][0] > 0.25]
    m.sort(key = lambda x: np.abs(x[1][1] - 0.5))
    if len(m) == 0:
        i = np.argmax(weights)
        m.append((i, means[i]))
    S = m[0][0]
    avg_RDR = sum([RDR[i] for i in range(len(clust)) if clust[i] == S]) / np.count_nonzero(clust == S)#clust.count(S)
    candidates = [(2**(i+1))/avg_RDR for i in range(max_wgd + 1)]
    
    return candidates, S

def estimate_ploidy_RC(RDR, states):
    X = RDR.to_numpy()

    candidates = []
    for s in states:
        est = []
        c = s - 0.1
        while c <= s + 0.1:
            est.append(np.sum(np.square(X*c - np.round(X*c))))
            c += 0.01
        val = min(est)
        candidates.append(val)
    return candidates

def gen_CN_pairs(maxCN):
    return [i for i in combinations_with_replacement([i for i in range(maxCN+1)], 2) if sum(i) <= maxCN]

# cell-specific, cluster-specific
def find_best_CN(RDR, BAF, scale, cov, pairs):
    size = len(RDR)
    best_pair, best_llh = None, None
    for CN1, CN2 in pairs:
        x = (CN1 + CN2) / scale
        if CN1 == 0 and CN2 == 0:
            y = 0
        else:
            y = min(CN1, CN2) / (CN1 + CN2)
        #print(cov)
        D = multivariate_normal(mean=[x, y], cov=cov)
        llh = sum([np.log(max(D.pdf([RDR[i], BAF[i]]), 1e-10)) for i in range(size)])
        if best_llh == None or llh > best_llh:
            best_pair = (CN1, CN2)
            best_llh = llh
    return best_pair, best_llh

#def check_zero_clust(mean):

# cell-specific
def compute_likelihood(clust, RDR, BAF, scale, n_wgd, S, pairs, covs):
    RDR_d, BAF_d = {}, {}
    for i in range(len(clust)):
        if clust[i] not in RDR_d:
            RDR_d[clust[i]] = []
            BAF_d[clust[i]] = []
        RDR_d[clust[i]].append(RDR[i])
        BAF_d[clust[i]].append(BAF[i])
    
    allele_CNs = {}
    total_llh = 0

    s_pair = (1 + n_wgd, 1 + n_wgd)
    alt_pairs = [p for p in pairs]

    best_pair, llh = find_best_CN(RDR_d[S], BAF_d[S], scale, covs[S], [s_pair])
    allele_CNs[S] = best_pair
    total_llh += llh
    del RDR_d[S]

    for c in RDR_d.keys():
        best_pair, llh = find_best_CN(RDR_d[c], BAF_d[c], scale, covs[c], alt_pairs)
        allele_CNs[c] = best_pair
        total_llh += llh
    
    return allele_CNs, total_llh



def get_cluster_CNs(cellnames, clust, RDR, BAF, means, covs, max_wgd, maxCN_user, ploidy_RDR):
    cell_clust = dict(zip(cellnames, np.split(clust, len(cellnames))))
    best_CNs = {}
    chosen = {}
    for cell in cellnames:
        cur_clust, cur_RDR, cur_BAF = cell_clust[cell], RDR.xs(cell), BAF.xs(cell)
        cur_weights = np.bincount(cur_clust) / len(cur_clust)
        while len(cur_weights) != len(means):
            cur_weights = np.append(cur_weights, 0)

        candidates, S = gen_candidates(cur_clust, means, cur_weights, max_wgd, cur_RDR)

        if ploidy_RDR:
            avg_RDR = sum([cur_RDR[i] for i in range(len(cur_clust)) if cur_clust[i] == S]) / np.count_nonzero(cur_clust == S)

            w = estimate_ploidy_RC(cur_RDR, candidates)
            n_wgd = np.argmin(w)
            scale = (2**(n_wgd+1))/avg_RDR
            maxCN_poss = np.ceil(scale*max(cur_RDR))
            maxCN = int(min(maxCN_poss, maxCN_user)) + 1
            pairs = gen_CN_pairs(maxCN)
            cur_CNs, cur_llh = compute_likelihood(cur_clust, cur_RDR, cur_BAF, scale, n_wgd, S, pairs, covs)
            chosen[cell] = n_wgd
            best_CNs[cell] = cur_CNs

        else:
            w = estimate_ploidy_RC(cur_RDR, candidates)
            num_CNs, candidate_CNs, llhs, weights = [], [], [], []
            for n_wgd, scale in enumerate(candidates):
                maxCN_poss = np.ceil(scale*max(cur_RDR))
                maxCN = int(min(maxCN_poss, maxCN_user)) + 1
                pairs = gen_CN_pairs(maxCN)
                num_CNs.append(len(pairs))
                cur_CNs, cur_llh = compute_likelihood(cur_clust, cur_RDR, cur_BAF, scale, n_wgd, S, pairs, covs)
                candidate_CNs.append(cur_CNs)
                llhs.append(cur_llh)
            m = len(cur_clust)
            best = np.argmin([w[i]*np.log(m)*num_CNs[i] - 2*llhs[i] for i in range(len(candidates))])
            chosen[cell] = best
            best_CNs[cell] = candidate_CNs[best]
            
    return best_CNs, chosen


def allele_caller(cell_names, best_CNs, clust):
    cell_clust = dict(zip(cell_names, np.split(clust, len(cell_names))))
    allele_CNs = {}
    for cell in cell_names:
        allele_CNs[cell] = [best_CNs[cell][i] for i in cell_clust[cell]]
    return allele_CNs

def get_segs(cell_names, clust, ensemble_bkpts, RDR, BAF, means, num_bin_per_chrom, min_seg_length):
    cell_clust = dict(zip(cell_names, np.split(clust, len(cell_names))))
    num_bins = len(cell_clust[cell_names[0]])
    allele_CNs = {}

    chrom_boundaries = [num_bin_per_chrom[0]]
    for x in num_bin_per_chrom[1:]:
        chrom_boundaries.append(x + chrom_boundaries[-1])

    ensemble_segments = {cell: np.split(cell_clust[cell], ensemble_bkpts[cell]) for cell in cell_names}
    distances = {}

    for cell in cell_names:
        cur_RDR, cur_BAF = RDR.xs(cell).to_numpy(), BAF.xs(cell).to_numpy()
        clusters_present = set(cell_clust[cell])
        seg_data = np.split(np.column_stack((cur_RDR, cur_BAF)), ensemble_bkpts[cell])

        def check_seg_sizes(seg_data):
            idxs = [i for i in range(len(seg_data))]
            idxs.sort(key=lambda x: len(seg_data[x]))
            best_idx, l = idxs[0], len(seg_data[idxs[0]])
            bin_start = 0
            for i in range(best_idx+1):
                bin_start += len(seg_data[i])
            return idxs[0], len(seg_data[idxs[0]])
        
        def get_bin_start(seg_data, idx):
            bin_start = 0
            for i in range(idx+1):
                bin_start += len(seg_data[i])
            return bin_start
        
        if min_seg_length > 1:
            idx, l = check_seg_sizes(seg_data)
            while l < min_seg_length:
                if idx == 0:
                    ensemble_bkpts[cell] = np.delete(ensemble_bkpts[cell], 0)
                elif idx == len(ensemble_bkpts[cell]):
                    ensemble_bkpts[cell] = np.delete(ensemble_bkpts[cell], -1)
                else:
                    start = idx - 1
                    end = idx
                    bin_s = get_bin_start(seg_data, start)
                    bin_e = get_bin_start(seg_data, end)
                    if bin_s in chrom_boundaries:
                        ensemble_bkpts[cell] = np.delete(ensemble_bkpts[cell], end)
                    elif bin_e in chrom_boundaries:
                        ensemble_bkpts[cell] = np.delete(ensemble_bkpts[cell], start)
                    else:

                        L = np.append(seg_data[idx - 1], seg_data[idx], axis=0)
                        R = np.append(seg_data[idx + 1], seg_data[idx], axis=0)
                        minL, minR = 1e6, 1e6
                        for j,m in enumerate(means):
                            if j in clusters_present:
                                x0,y0 = euc(np.mean(L, axis=0), m), euc(np.mean(R, axis=0), m)
                                if x0 < minL:
                                    minL = x0
                                if y0 < minR:
                                    minR = y0
                        if minL <= minR:
                            ensemble_bkpts[cell] = np.delete(ensemble_bkpts[cell], start)
                        else:
                            ensemble_bkpts[cell] = np.delete(ensemble_bkpts[cell], end)
                        
                seg_data = np.split(np.column_stack((cur_RDR, cur_BAF)), ensemble_bkpts[cell])
                idx, l = check_seg_sizes(seg_data)
            
            ensemble_segments[cell] = np.split(cell_clust[cell], ensemble_bkpts[cell])
    
        distances[cell] = [[] for i in range(len(ensemble_segments[cell]))]
        seg_data = np.split(np.column_stack((cur_RDR, cur_BAF)), ensemble_bkpts[cell])
        for i,seg in enumerate(seg_data):
            seg_mean = np.mean(seg, axis=0)
            for j,m in enumerate(means):
                if j in clusters_present:
                    eucd = euc(seg_mean, m)
                    distances[cell][i].append((j, eucd))

    cell_assignments = {cell: [] for cell in cell_names}
    for cell in cell_names:
        for i, dists in enumerate(distances[cell]):
            dists.sort(key = lambda x: x[1])
            k, eucd = dists.pop(0)
            seg = ensemble_segments[cell][i]
            cell_assignments[cell].extend([k for j in range(len(seg))])
    
    new_clusts = []
    for cell in cell_names:
        new_clusts.extend(cell_assignments[cell])

    return np.array(new_clusts)