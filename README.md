# SEACON
SEACON (Single-cell Estimation of Allele-specific COpy Numbers) is a tool for allele-specific copy number profiling from single-cell DNA-sequencing data.

SEACON can be cited as follows:
"Improved allele-specific single-cell copy number estimation in low-coverage DNA-sequencing".
Samson Weiner, Bingjun Li, and Sheida Nabavi.
Under Review.
 
## Installation
SEACON is written in python and R and requires interpreters for both languages. To install, first clone the repository:
```
git clone https://github.com/NabaviLab/SEACON.git
```
Next, finish setting up your environment by installing the required packages and downloading preprocessing files.

### Dependencies
SEACON requires the following python packages be installed:
* [Numpy](https://numpy.org/)
* [Pandas](https://pandas.pydata.org/)
* [Scipy](https://scipy.org/)
* [scikit-learn](https://scikit-learn.org/stable/)
* [rpy2](https://rpy2.github.io/)
* [Pysam](https://pysam.readthedocs.io/en/latest/api.html)
* [Pyfaidx](https://github.com/mdshw5/pyfaidx)

The following R packages are also required:
* [dnacopy](https://bioconductor.org/packages/release/bioc/html/DNAcopy.html)

Additionally, SEACON requires the following bioinformatics tools be available:
* [bedtools](https://bedtools.readthedocs.io/en/latest/)
* [bcftools](https://samtools.github.io/bcftools/bcftools.html)
* [bigWigAverageOverBed](https://www.encodeproject.org/software/bigwigaverageoverbed/)

### Preparation
For data preprocessing, SEACON requires a few additional files. First, a reference genome must be provided. This will typically be hg38 or hg19, and can be obtained from the UCSC genome browser.

Second, a mappability file in bigWig format must be provided for the given reference genome. These can be obtained from the same organization (for example, [here](https://genome.ucsc.edu/cgi-bin/hgFileUi?db=mm9&g=wgEncodeMapability)). The authors of the SeCNV have also made these files for hg38 and hg19 easily accessable on their [Google Drive](https://drive.google.com/drive/folders/1XGuXa9muRtOiAjtz1J4gWO5Qk3c5-PwS). 

SEACON also requires precomputed BAFs. For this, we use [CHISEL](https://github.com/raphael-group/chisel). We recommend installing CHISEL in a fresh conda environment following the instructions at the link provided. We have provided a script in the `scripts` folder for generating the BAFs using CHISEL.

## Usage
SEACON is run by calling the `seacon.py` script in the *core* directory with a selected mode as follows:
```
python core/seacon.py [mode] [-r reference] [-b bin size]...
```
There are a total 5 modes: *prep*, *norm*, *pseudonormal*, *call*, and *full*. Using the *full* mode runs the full pipeline of SEACON assuming all the inputs are prepared. However, computing BAFs with CHISEL or other tools may require SNPs be called from a pseudonormal sample. In that case, SEACON can be used to first build a pseudonormal bam file, and then the rest of the pipeline can be completed later.

The following section details how to run SEACON using each mode.

**prep mode.** The *prep* mode of SEACON partitions the genomeinto bins, then extracts the raw readcounts, mappability, and GC content values for each bin. It is executed as follows:
```
python core/seacon.py prep -i bam_dir -o out_dir -r hg38.fa -b 1000000 --map-file hg38_mappability.bigWig
```

**norm mode.** The *norm* mode of SEACON identifies normal cells in the sample, then normalizes the read counts. It is executed as follows:
```
python core/seacon.py norm -o out_dir --min-norm 10
```

**pseudnormal mode.** The *pseudnormal* mode of SEACON builds a BAM file from the normal cells identified in the previous step, and is executed as follows:
```
python core/seacon.py pseudnormal -i bam_dir -o out_dir
```
See the [subsection](###-Obtaining-list-of-phased-SNPs) under Miscelaneous to obtain the phased SNPs from the pseudonormal sample.

**call mode.** The *call* mode of SEACON is the segmentation and copy number calling part of the pipeline, and requires having run SEACON in *prep* and *norm* mode and the input BAFs be provided. It is executed as follows:
```
python core/seacon.py call -o out_dir --baf-path [path_to_bafs] --tolerance 0.2 --upper-filter 5 --max-wgd 2
```
The path to the input BAFs is required in this step. If the BAFs are generated with CHISEL, simply set *--baf-path* to the location of the directory created when running CHISEL. Otherwise, set to a tsv file containing the BAFs. The tsv file should have 5 columns, which in order should be: chromosome name, start index, end index, cell name, BAF value. There should also be a header row containing entries 'chrom', 'start', 'end', 'cell', and 'BAF'.

For a full list of parameters and their default, use the following command:
```
python core/seacon.py --help
```

## Miscelaneous

### Preprocessing
If raw reads are acquired in fastq format, we recommend the following steps (note that the bioinformatics tools *bwa* and *samtools* are required):
```
bwa mem -t 8 hg38.fa cell.fastq.gz > cell.sam
samtools view -q 40 -hbS -o cell.bam cell.sam
samtools sort -o cell.sorted.bam cell.bam
samtools index cell.sorted.bam
```

### Obtaining list of phased SNPs
If a matched-normal sample is available, or if a pseudonormal sample is constructed using the *pseudonormal* mode of SEACON, it is necessary to obtain a set of phased SNPs to compute the BAFs. For computing BAFs, we recommend using [CHISEL](https://github.com/raphael-group/chisel). Among other inputs, CHISEL requires a vcf file for the SNPs or optionally tsv file with 3 columnns: the first column is the chromosome name, the second column is the position, and the third is 0|1 or 1|0, indicating the phase of the SNP.

We first describe how to obtain this file. First, call SNPs using bcftools:
```
bcftools mpileup --skip-indels --ignore-RG -f hg38.fa pseudonormal.bam | bcftools call -mv -Ob -o snps.bcf
```

Next, phase the SNPs in the bcf file using Eagle2. This can be done using the [Michigan Imputation Server](), or users can run Eagle2 [locally](https://alkesgroup.broadinstitute.org/Eagle/) if a reference panel is downloaded locally, for example from [this reference](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/) from 1000 Genomes. 

For running Eagle2 through the Michigan Imputation Server, follow the directions given in the documentation. If using the tool locally, run on each chromosome as follows:
```
for chr in {1..22}; do
    mkdir snps phased_snps
    bcftools view -Ob -r chr${chr} -o snps/chr${chr}.bcf.gz snps.bcf
    bcftools index snps/chr${chr}.bcf.gz
    ./Eagle2/eagle --vcfTarget snps/chr${chr}.bcf --vcfRef panel/chr${chr}.bcf.gz --geneticMapFile \\
    Eagle2/tables/genetic_map_hg38_withX.txt.gz --outPrefix phased_snps/chr${chr} --vcfOutFormat b
done
```

Next, concatenate all chromosomes into a single file.
```
bcftools concat -Ob -o all_phased.bcf phased_snps/chr{1..22}.bcf
```
Filter for heterozygous sites:
```
bcftools view -Ob -g het -o all_phased_het.bcf all_phased.bcf
bcftools view -Oz -o all_phased_het.vcf.gz all_phased_het.bcf
```

Finally, use the `write_snps.py` script located in the *scripts* directory to obtain a tsv file as follows (note the python package vcf is required):
```
python write_snps.py -i all_phased_het.vcf.gz -o phases.tsv
```

### Generating BAFs with CHISEL
Using the *phases.tsv* file generated in the previous step, users can run CHISEL to generate the input BAFs using the `run_chisel.sh` script located in the *scripts* directory. Ensure the chisel environment is active when running the script.

When running SEACON, 

## Simulations
Datasets were generated using the [CNAsim](https://github.com/samsonweiner/CNAsim) simulator. Information for recreating datasets used in the paper can be found in the README.md file in the *simulation* directory.