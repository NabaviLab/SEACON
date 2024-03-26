import vcf
import argparse

def write_phase_single(inp_vcf, out_file):
    inp = vcf.Reader(filename=inp_vcf)
    with open(out_file, 'w+') as f:
        for rec in inp:
            chrom, pos = rec.CHROM, rec.POS
            GT = rec.genotype('local0/pseudonormal.bam').data.GT
            f.write(f'{chrom}\t{pos}\t{GT}\n')

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--in-path', type=str, default=None, help='Input vcf.gz file with phase heterozygous SNPs.')
parser.add_argument('-o', '--out-path', type=str, default=None, help='Output path.')
args = parser.parse_args()

write_phase_single(args.in_path, args.out_path)