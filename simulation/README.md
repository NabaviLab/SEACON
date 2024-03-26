# Simulation README

## Direct read count/BAF simulations

To generate datasets with read and BAFs generated directly, first generate the copy number profiles using CNAsim, then run the `gen_direct_data.py` script.

### Command to run CNAsim for generating direct read count/BAF data:

```
python CNAsim/main.py -m 0 -o seg_data/d1/sim_output -n 100 -c 3 -c2 20 -N 10 -v -i3 1.5
python CNAsim/main.py -m 0 -o seg_data/d2/sim_output -n 100 -c 3 -c2 20 -N 10 -v -i3 1.5 -l 10000000 -l1 5000000
python CNAsim/main.py -m 0 -o seg_data/d3/sim_output -n 100 -c 3 -c2 20 -N 10 -v -i3 1.5 -B 250000 -l 10000000 -l1 5000000
```


## Synthetic sequence read simulations

To use CNAsim to generate synthetic sequencing reads datasets, run in mode 2 with the necessary parameters. An example run for dataset A2 is shown below:

```
python3 CNAsim/main.py -o A2 -m 2 -r1 ../alt/genomes/alt1.fa -r2 ../alt/genomes/alt2.fa -n 100 -n1 0.25 -c 3 -c2 20 -U -v -q 0.5 -u 0.5 -i 3 -i1 3 -i2 3 -C 0.1 -M -B 1000000 -W 1000000 -k 10000 -l1 10000
```

Then, apply standard alignment/quality control steps to obtain BAM files (details are in the preprocessing section of main README file).

If evaluating allele frequencies, it is imperative to use two haploid reference genomes (in the command called alt1.fa and alt2.fa) where each has a unique set of SNPs. Details on how these references were constructed is given in the main SEACON paper.
