# Generating CNA profiles
Details of how the CNA profiles were generated for each dataset are provided below.

## T10 and T16
Raw read counts in the form of fastq files were obtained from the Sequence Read Archive (SRA) under accession code SRA018951.105. We then ran bwa alignment, filtered out low quality reads, and indexed the resulting bam files following the raw data preprocessing steps provided in the README file at the root of the repository.

SCOPE profiles were generated following the instructions provided in the user manual with a bin length of 1Mbp. 

To generate the SEACON profiles, we first ran  `seacon prep_readcount` followed by `seacon pseudonormal` to identify normal cells and create a pseudonormal bam file. Next, we obtained a set of phased SNPs from the pseudonormal file following the instructions provided in the README file at the root of the repository. Next, we called the BAF estimation routine of CHISEL and `seacon prep_baf` to prepare the input BAFs. Finally, we ran the main call function of seacon using the following command:
```
seacon call -o SEACON_T16 --upper-filter 5 --tolerance 0.15 --max-wgd 1
```

## 10x
Processed CNA profiles of CHISEL as well as input RDRs and BAFs for all cells were acquired from the [chisel data repository](https://github.com/raphael-group/chisel-data). We then filtered the input RDRs and BAFs for cells from section E and fed these directly into SEACON. The main call function was executed as follows:
```
seacon call -o SEACON_10x_E --upper-filter 20 --tolerance 0.05
```