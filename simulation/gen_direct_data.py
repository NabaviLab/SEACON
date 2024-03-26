import numpy as np
import pandas as pd
import os

def read_profiles(dir_path):
    out_dir = os.path.join(dir_path, 'profiles.tsv')
    df = pd.read_csv(out_dir, delimiter='\t')
    df = df.set_index(['CELL', 'chrom', 'start'])
    df = df.drop(columns=['end'])
    df = df.sort_index()
    df = df['CN states'].str.split(',', expand=True)
    return df

def collect_breakpoints(df, max_CN=6):
    cell_names = list(df.index.get_level_values(0).drop_duplicates())
    chrom_names = list(df.index.get_level_values(1).drop_duplicates())
    cell_names.sort(key = lambda x: int(x[4:]))
    chrom_names.sort(key = lambda x: int(x[3:]))

    bin_to_id = {chrom: {} for chrom in chrom_names}
    id_to_bin = {}
    counter = 0
    for chrom in chrom_names:
        for start in df.xs(chrom, level='chrom').index.get_level_values(1).drop_duplicates():
            bin_to_id[chrom][start] = counter
            id_to_bin[counter] = (chrom, start)
            counter += 1
    
    CNPs, breakpoints = {}, {}
    for cell in cell_names:
        CNPs[cell] = [[], []]
        cur_df = df.xs(cell, level='CELL')
        bin_ids = list(id_to_bin.keys())
        bin_ids.sort()
        for b in bin_ids:
            chrom, start = id_to_bin[b]
            CNs = cur_df.loc[chrom, start]
            CNPs[cell][0].append(min(int(CNs[0]), max_CN))
            CNPs[cell][1].append(min(int(CNs[1]), max_CN))

        breakpoints[cell] = []
        prev_a, prev_b = CNPs[cell][0][0], CNPs[cell][1][0]
        for b in bin_ids[1:]:
            cur_a, cur_b = CNPs[cell][0][b], CNPs[cell][1][b]
            if cur_a != prev_a or cur_b != prev_b:
                breakpoints[cell].append(b)
                prev_a, prev_b = cur_a, cur_b
    
    return cell_names, id_to_bin, CNPs, breakpoints

def gen_readcount_data(cell_names, bin_ids, CNPs, scale=100):
    readcounts, bafs = {}, {}
    for cell in cell_names:
        readcounts[cell], bafs[cell] = {'low': [], 'med': [], 'high': []}, {'low': [], 'med': [], 'high': []}

        for b in bin_ids:
            a, b = CNPs[cell][0][b], CNPs[cell][1][b]
            tot_CN = a + b
            if tot_CN == 0:
                exp_rc = np.sqrt(scale)*1.2
                rat = 0.5
            else:
                exp_rc = tot_CN*scale
                rat = np.random.choice([a,b]) / tot_CN

            low_rc = np.random.poisson(exp_rc)
            med_rc = round(low_rc * np.random.normal(1, 0.1))
            high_rc  = round(low_rc * np.random.normal(1, 0.25))

            low_baf = np.random.normal(rat, 0.01)
            med_baf = np.random.normal(rat, 0.02)
            high_baf = np.random.normal(rat, 0.03)

            low_baf = np.abs(low_baf)
            if low_baf > 1:
                low_baf = 2 - low_baf
            med_baf = np.abs(med_baf)
            if med_baf > 1:
                med_baf = 2 - med_baf
            high_baf = np.abs(high_baf)
            if high_baf > 1:
                high_baf = 2 - high_baf

            readcounts[cell]['low'].append(low_rc)
            readcounts[cell]['med'].append(med_rc)
            readcounts[cell]['high'].append(high_rc)

            bafs[cell]['low'].append(low_baf)
            bafs[cell]['med'].append(med_baf)
            bafs[cell]['high'].append(high_baf)
    
    return readcounts, bafs

def create_dataset(sim_path, scale, skip_gt=False):
    scale_path = os.path.join(sim_path, 'scale' + str(scale))
    if not os.path.isdir(scale_path):
        os.makedirs(scale_path)
    print('Collecting CNPs and computing breakpoints')
    df = read_profiles(os.path.join(sim_path, 'sim_output'))
    cell_names, id_to_bin, CNPs, breakpoints = collect_breakpoints(df)
    bin_ids = list(id_to_bin.keys())
    bin_ids.sort()
    print('Generating synthetic read counts and bafs')
    readcounts, bafs = gen_readcount_data(cell_names, bin_ids, CNPs, scale=scale)

    print('Writing data')
    headers = ['cell'] + [str(b) for b in bin_ids]
    #levels = ['low', 'med', 'high']
    levels = ['low', 'med']

    if not skip_gt:
        f = open(os.path.join(sim_path, 'CNP.tsv'), 'w+')
        f.write('\t'.join(headers) + '\n')
        for cell in cell_names:
            f.write('\t'.join([cell] + [str(CNPs[cell][0][b]) + ',' + str(CNPs[cell][1][b]) for b in bin_ids]) + '\n')
        f.close()

    for level in levels:
        f = open(os.path.join(scale_path, 'readcounts_' + level + '.tsv'), 'w+')
        f.write('\t'.join(headers) + '\n')
        for cell in cell_names:
            f.write('\t'.join([cell] + [str(readcounts[cell][level][b]) for b in bin_ids]) + '\n')
        f.close()

        f = open(os.path.join(scale_path, 'bafs_' + level + '.tsv'), 'w+')
        f.write('\t'.join(headers) + '\n')
        for cell in cell_names:
            f.write('\t'.join([cell] + [str(bafs[cell][level][b]) for b in bin_ids]) + '\n')
        f.close()

    if not skip_gt:
        f = open(os.path.join(sim_path, 'breakpoints.tsv'), 'w+')
        for cell in cell_names:
            f.write('\t'.join([cell] + [str(i) for i in breakpoints[cell]]) + '\n')
        f.close()


# First parameter is path to CNAsim directory, second parameter is where to put data
#print('d6')
create_dataset('d1', 100)
create_dataset('d1', 300, skip_gt=True)
create_dataset('d1', 600, skip_gt=True)

create_dataset('d2', 100)
create_dataset('d2', 300, skip_gt=True)
create_dataset('d2', 600, skip_gt=True)

create_dataset('d3', 100)
create_dataset('d3', 300, skip_gt=True)
create_dataset('d3', 600, skip_gt=True)