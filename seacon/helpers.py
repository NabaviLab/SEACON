import numpy as np
import pandas as pd
import os

# Assumes data has been processed with read_data from complete_format
def flatten(df, cell_names):
    chrom_names = list(df.index.get_level_values(1).drop_duplicates())
    #chrom_names.sort(key = lambda x: int(x[3:]))
    flat_data = {}
    for cell in cell_names:
        cell_df = df.xs(cell, level='cell')
        flat_data[cell] = [i for chrom in chrom_names for i in list(cell_df.xs(chrom, level='chrom').unstack(level='start').to_numpy())]
    flat_df = pd.DataFrame.from_dict(flat_data, orient='index')
    flat_df.index.name = 'cell'
    return flat_df

def unflatten(df, bin_ids, col_name='count'):
    cell_names = list(df.index)
    frames = []
    for cell in cell_names:
        cell_df = df.xs(cell)
        cell_frames = [[cell, x[0], x[1], x[2], cell_df[i]] for i,x in enumerate(bin_ids)]
        frames.append(pd.DataFrame(cell_frames, columns=['cell', 'chrom', 'start', 'end', col_name]))
    df = pd.concat(frames)
    df = df.set_index(['cell', 'chrom', 'start'])
    return df

def get_coords(df):
    chrom_names = list(df.index.get_level_values(1).drop_duplicates())
    chrom_names.sort(key = lambda x: int(x[3:]))
    coords = {chrom: list(df.xs(chrom, level='chrom').index.get_level_values(1).drop_duplicates()) for chrom in chrom_names}
    return coords

def get_bins_per_chrom(coords):
    num_bins_per_chrom = []
    prev_chrom = None
    for x in coords:
        chrom = x[0]
        if chrom != prev_chrom:
            prev_chrom = chrom
            num_bins_per_chrom.append(0)
        num_bins_per_chrom[-1] += 1
    return num_bins_per_chrom

def get_chrom_start_index(coords, chrom_names):
    chrom_names.sort(key = lambda x: int(x[3:]))
    chrom_start_bins = []
    idx = 0
    for chrom in chrom_names[:-1]:
        idx += len(coords[chrom])
        chrom_start_bins.append(idx)
    return chrom_start_bins

# Maps bins across all chroms to integer ids from region file or coords df.
def index_bins(region_path, coords = None):
    cur_id = 0
    bin_ids = []
    if not coords:
        with open(region_path) as f:
            for line in f:
                info = line.strip().split('\t')
                chrom, start, end = info[0], int(info[1]), int(info[2])
                bin_ids.append((chrom, start, end))
                #cur_id += 1
    else:
        ends = {}
        with open(region_path) as f:
            for line in f:
                info = line.strip().split('\t')
                chrom, start, end = info[0], int(info[1]), int(info[2])
                if chrom not in ends:
                    ends[chrom] = {}
                ends[chrom][start] = end
        chrom_names = list(coords.keys())
        chrom_names.sort(key = lambda x: int(x[3:]))
        for chrom in chrom_names:
            for start in coords[chrom]:
                bin_ids.append((chrom, start, ends[chrom][start], cur_id))
                cur_id += 1
    return bin_ids

def save_flat_data(df, outdir, name, bin_ids = None):
    if bin_ids:
        full_df = unflatten(df, bin_ids)
        full_df.to_csv(os.path.join(outdir, name), sep='\t')
    else:
        df.to_csv(os.path.join(outdir, name), sep='\t')
    
def combine_flat_data(read_counts_df, baf_df):
    binID_to_cell = {}
    cell_bin_to_binID = {}
    cell_names = list(read_counts_df.index.get_level_values(0).drop_duplicates())
    num_bins = len(read_counts_df.columns)
    r = read_counts_df.to_numpy().flatten()
    b = baf_df.to_numpy().flatten()
    #b = np.minimum(b, 1-b)
    #f(b)
    flat_data = list(zip(r,b))
    for i, cell in enumerate(cell_names):
        cell_bin_to_binID[cell] = []
        for j in range(num_bins):
            binID_to_cell[i*num_bins + j] = [cell, j]
            cell_bin_to_binID[cell].append(i*num_bins + j)
    return flat_data, binID_to_cell, cell_bin_to_binID

# returns dictionary of cells with breakpoints
def compute_bkpts(cell_names, binID_to_cell, cluster_assignments):
    bkpts = {}
    for cell in cell_names:
        bkpts[cell] = []
    cur_cell, cur_clust = binID_to_cell[0][0], cluster_assignments[0]
    for i,c in enumerate(cluster_assignments[1:]):
        if binID_to_cell[i][0] != cur_cell:
            cur_cell = binID_to_cell[i][0]
            cur_clust = c
        else:
            if c != cur_clust:
                bkpts[cur_cell].append(binID_to_cell[i][1]+1)
                cur_clust = c
    return bkpts

def read_breakpoints(filename):
    breakpoints = {}
    with open(filename) as f:
        for line in f.readlines():
            cell, *bkpts = line.strip().split('\t')
            breakpoints[cell] = [int(i) for i in bkpts]
    return breakpoints

def save_breakpoints(bkpts, cell_names, filename):
    f = open(filename, 'w+')
    for cell in cell_names:
        f.write('\t'.join([cell] + [str(i) for i in bkpts[cell]]) + '\n')
    f.close()

def save_clusters(cell_names, binID_to_cell, cluster_assignments, filename):
    clusters = {}
    for cell in cell_names:
        clusters[cell] = []
    for i,c in enumerate(cluster_assignments):
        clusters[binID_to_cell[i][0]].append(c)
    f = open(filename, 'w+')
    for cell in cell_names:
        f.write('\t'.join([cell] + [str(i) for i in clusters[cell]]) + '\n')
    f.close()

def get_RDR_CHISEL(combo_path, cellmap = None):
    df = pd.read_csv(combo_path, sep='\t', names=['chrom', 'start', 'end', 'barcode', 'normcount', 'readcount', 'RDR', 'Acount', 'Bcount', 'BAF'])
    df = df.drop(['normcount', 'Acount', 'Bcount', 'BAF'], axis=1)
    df = df.set_index(['barcode', 'chrom', 'start'])

    if cellmap:
        df = df.rename(index=cellmap)
    df.index.names = ['cell', 'chrom', 'start']

    #readcount_df = df.drop(['RDR'], axis=1)
    #RDR_df = df.drop(['readcount'], axis=1)
    #return readcount_df, RDR_df
    return df

# stored in the 'locations.tsv' file
def get_coords_SCOPE(location_path):
    loc_df = pd.read_csv(location_path)
    loc_df = loc_df.drop(['end', 'width', 'strand', 'gc', 'mapp'], axis=1)
    loc_df['start'] = loc_df['start'].astype(int)
    loc_df['start'] = loc_df['start'] - 1

    chrom_names = list(loc_df['seqnames'].drop_duplicates())
    coords = {chrom: list(loc_df.loc[loc_df['seqnames'] == chrom]['start']) for chrom in chrom_names}
    return coords

# stored in 'read.csv' file
def get_RDR_SCOPE(path):
    f = open(os.path.join(path, 'sample.txt'))
    samples = [i[:-1] for i in f.readlines()]
    f.close()

    locations = pd.read_csv(os.path.join(path, 'locations.csv'), delimiter=',')
    output = pd.read_csv(os.path.join(path, 'read.csv'))
    output.columns = samples

    avgs = {}
    for cell in samples:
        avgs[cell] = np.mean(output[cell])

    frames = []
    for i in locations.index:
        loc_info = locations.loc[i,]
        chrom, start, end = loc_info[0], int(loc_info[1])-1, loc_info[2]
        RCs = output.loc[i,]
        for j, cell in enumerate(samples):
            rc = RCs[j]
            row = [cell, chrom, start, end, rc, rc/avgs[cell]]
            frames.append(row)
    df = pd.DataFrame(frames, columns=['cell', 'chrom', 'start', 'end', 'readcount', 'RDR'])
    df = df.set_index(['cell', 'chrom', 'start'])
    return df

# stored in 'genome_cov.bed' file
def get_RDR_SeCNV(readcount_path):
    df = pd.read_csv(readcount_path, sep='\t')
    cellnames = list(df.columns)[3:-1]

    df = pd.melt(df, id_vars=['chromosome', 'start', 'stop'], value_vars=cellnames, var_name='cell', value_name='readcount')
    df['start'] = df['start'].astype(int)
    df['start'] = df['start'] - 1
    df = df.set_index(['cell', 'chromosome', 'start'])
    df.index.names = ['cell', 'chrom', 'start']
    df = df.rename(columns={'stop': 'end'})

    avgs = {}
    for cell in cellnames:
        cell_df = df.xs(cell, level='cell')
        avg_rc = np.mean(cell_df['readcount'])
        avgs[cell] = avg_rc
    df['avgs'] = df.apply(lambda x: avgs[x.name[0]], axis=1)
    df['RDR'] = df['readcount'] / df['avgs']
    df = df.drop(['avgs'], axis=1)
    return df

# stored in the 'ref.bed' file
def get_coords_seCNV(location_path):
    coords = {}
    with open(location_path) as f:
        for line in f:
            info = line.strip().split('\t')
            chrom, start = info[0], int(info[1])-1
            if chrom not in coords:
                coords[chrom] = []
            coords[chrom].append(start)
    return coords


def correct_RC_BAF(RC_df, BAF_df, bins1, bins2, normal_cells=[]):
    cell_names1, cell_names2 = list(RC_df.index), list(BAF_df.index)
    cell_names = list(set(cell_names1).intersection(cell_names2))
    if len(normal_cells) > 0:
        for c in normal_cells:
            if c not in cell_names2:
                cell_names.append(c)
                BAF_df.loc[c] = [0.5 for i in range(len(bins2))]
        
    RC_df, BAF_df = RC_df[RC_df.index.isin(cell_names)], BAF_df[BAF_df.index.isin(cell_names)]

    #shared = list(set(bins1).intersection(set(bins2)))
    shared = [b for b in bins1 if b in bins2]
    keep1 = [i for i in range(len(bins1)) if bins1[i] in bins2]
    keep2 = [i for i in range(len(bins2)) if bins2[i] in bins1]
    RC_df_filt = RC_df.filter(items=keep1, axis=1)
    BAF_df_filt = BAF_df.filter(items=keep2, axis=1)

    RC_df_filt = RC_df_filt.rename(columns=dict(zip(RC_df_filt.columns, range(len(shared)))))
    BAF_df_filt = BAF_df_filt.rename(columns=dict(zip(BAF_df_filt.columns, range(len(shared)))))

    return RC_df_filt, BAF_df_filt, shared, cell_names