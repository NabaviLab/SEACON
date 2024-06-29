# Scripts instructions

## Instructions for using the breakpoint filtering heuristic with CHISEL
The `smooth_chisel.py` script can be used to smooth the results of chisel by using the breakpoint filtering of CHISEL. Having previously run both SEACON and CHISEL, the script is used as follows:
```
python smooth_chisel.py -i path_to_seacon_dir -c path_to_chisel_dir -o path_to_output_dir --upper-filter 5
```

## Instructions for using breakpoint filtering heuristic with any output CNPs
The `smooth_general.py` script can be used to smooth the results of any output CNPs. The CNPs should be in tsv format, where each row is the cell name followed by the CNs across the genome. The header row should be numbered 0-n separated by tabs, where n is the number of bins. The breakpoint file should be in the same format, except with no header row and the values are the bin indexes where breakpoints occur (see the local_bkpts.tsv file in the SEACON output). Then the script is used as follows:
```
python smooth_general.py -i path_to_input_cns -b path_to_input_bkpts -o path_to_output_dir --upper-filter 5
```

## Instructions for getting phased SNPs from bam file
The `gen_vcf.sh` script can be used to call and phase SNPs from a normal or pseudo-normal sample in bam format. This requires bcftools and a phasing tool are installed. By default, it uses eagle2. Simply alter the paths at the top of the file to the correct locations, and run the script with
```
bash gen_vcf.sh
```