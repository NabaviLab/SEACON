import os

def check_prep_readcount(args):
    if args['reference'] is None:
        print('Must provide a reference.')
        return -1
    if not os.path.exists(args['reference']):
        print('Cannot find reference file')
        return -1
    if args['bin_size'] is None:
        print('Bin size required')
        return -1
    if not os.path.exists(args['map_file']):
        print('Must provide bigWig mappability file.')
        return -1
    if not os.path.isdir(args['bam_path']):
        print('Must specify path to directory containing bam files.')
        return -1
    return 0

def check_prep_baf(args):
    if not os.path.isdir(args['bam_path']):
        print('Must specify path to directory containing bam files.')
        return -1
    if args['precomputed_baf']:
        if not os.path.exists(args['precomputed_baf']):
            print('Cannot find input baf file.')
            return -1
    elif args['vcf']:
        if not os.path.exists(args['vcf']):
            print('Cannot find input vcf file.')
            return -1
    else:
        print('Must provide input vcf with phased snps or precomputed BAFs.')
        return -1
    return 0