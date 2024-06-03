# SEACON
SEACON (Single-cell Estimation of Allele-specific COpy Numbers) is a tool for allele-specific copy number profiling from single-cell DNA-sequencing data.

SEACON can be cited as follows:
"Improved allele-specific single-cell copy number estimation in low-coverage DNA-sequencing".
Samson Weiner, Bingjun Li, and Sheida Nabavi.
Under Review.

## Contents
- [Installation](#Installation)
    - [Standard Installation](#standard-installation)
    - [Custom Installation](#custom-installation)
- [Inputs and Preparation](#inputs-and-preparation)
    - [Sequencing Reads](#sequencing-reads)
    - [Reference Genome](#reference-genome)
    - [Mappability File](#mappability-file)
    - [Phasing Information](#phasing-information)
- [Usage](#usage)
    - [Readcount Prep](#prep_readcount-mode)
    - [BAF Prep](#prep_baf-mode)
    - [CNA calling](#call-mode)
    - [Constructing pseudonormal](#pseudonormal-mode)
    - [Output files](#output-files)
- [Utilities](#utilities)
    - [Raw Data Preprocessing](#raw-data-preprocessing)
    - [Obtaining Phased SNPs](#obtaining-phased-snps)
    - [Example commands](#example-commands)
- [Simulations](#simulations)

## Installation
To install SEACON, begin by cloning the repository, then cd into it:
```
git clone --recursive https://github.com/NabaviLab/SEACON.git
cd SEACON
```
Next, follow either the instructions for either the standard installation or custom installation.

### Standard Installation
To reproduce the environment needed to use SEACON, with all its dependencies and versions, create a new conda environment from the provided **environment.yml** file:
```
conda env create -f environment.yml
```

This will create a new conda environment named SEACON. Next, execute the installation:
```
python setup.py install
```

Lastly, be sure to activate the environment before using the tool:
```
conda activate SEACON
```

### Custom Installation
The environment can also be assembled manually. SEACON reliably works for Python version 3.10 and R version 4.1. Additionally, there are a number of python dependencies, one R package, and a few standard bioinformatics tools. All are available on conda, and are listed below.

#### Dependencies
SEACON requires the following python packages are installed:
* [Numpy](https://numpy.org/)
* [Pandas](https://pandas.pydata.org/)
* [Scipy](https://scipy.org/)
* [scikit-learn](https://scikit-learn.org/stable/)
* [rpy2](https://rpy2.github.io/)
* [Pysam](https://pysam.readthedocs.io/en/latest/api.html)
* [Pyfaidx](https://github.com/mdshw5/pyfaidx)
* [PyVCF](https://pyvcf.readthedocs.io/en/latest/)

The following R packages are also required:
* [dnacopy](https://bioconductor.org/packages/release/bioc/html/DNAcopy.html)

Additionally, SEACON requires the following bioinformatics tools be available:
* [bedtools](https://bedtools.readthedocs.io/en/latest/)
* [bcftools](https://samtools.github.io/bcftools/bcftools.html)
* [bigWigAverageOverBed](https://www.encodeproject.org/software/bigwigaverageoverbed/)
These tools must be added to your `PATH` variable.

Additionally, SEACON integrates [CHISEL](https://github.com/raphael-group/chisel). CHISEL is written in python2.7 and thus causes dependency conflicts, so it should be installed in a separate conda environment. SEACON assumes it is installed in an environment named "SEACON_chisel". Note that this step is automated with the provided **initialize.sh** script.

## Inputs and Preparation

### Sequencing reads
SEACON takes as input a path to a directory containing bam files, one for each cell. It assigns the name of the file to be the cell id. Information on obtaining bam files from raw fastq files is provided in the [Raw data preprocessing](#raw-data-preprocessing).

### Reference genome
A reference genome is needed for some of the steps in SEACON. For human reference genomes, this will typically be hg38 or hg19, and can be obtained from the UCSC genome browser [here](https://hgdownload.soe.ucsc.edu/downloads.html#human).

### Mappability file
A mappability file in bigWig format must be provided for the given reference genome. These can be obtained from the same organization (for example, [here](https://genome.ucsc.edu/cgi-bin/hgFileUi?db=mm9&g=wgEncodeMapability)). The 100-mers mappability track for both hg38 and hg19 are easily accessable on this [Google Drive](https://drive.google.com/drive/folders/1XGuXa9muRtOiAjtz1J4gWO5Qk3c5-PwS) (Credit to Wang Ruohan). 

### Phasing information
To compute BAFs, SEACON uses the estimation procedure of [CHISEL](https://github.com/raphael-group/chisel), and is included as a submodule in this repository. This requires an input vcf file containing phased SNPs from a matched-normal sample or a pseudo-normal sample built from identified normal cells (see [this section](#pseudonormal-mode)). Scripts are provided for obtaining the phased SNPs using additional bioinformatics tools. More information appears in the section [Obtaining phased SNPs](#obtaining-phased-snps) below. Additionally, precomputed BAFs can be provided if estimated using a different algorithm. For convenience, if using the BAF estimation procedure of CHISEL, the outputs of CHISEL will also be generated for comparison.


## Usage
SEACON is run by calling the command `seacon` with a selected mode as follows (or by running the *core/seacon.py* script):
```
seacon [mode] [-o outdir] [options]
```
The modes are as follows: `prep_readcount`, `prep_baf`, `call`, and `pseudonormal`. SEACON is intended to be run using the first three modes in that order, which can also be achieved in a single call by passing mode `full`. The `pseudonormal` mode is designed to build a pseudobulk bam file from the normal cells to facilitate the BAF computation, but is only necessary if a matched-normal sample is unavailable. Each of the four modes are explained in further detail below.

For a full list of parameters and their default, use the following command:
```
seacon --help
```
Additionally, please be advised that multiple cores can be utilized by specifying the `--num-processors` parameter.

### prep_readcount mode
The *prep_readcount* mode partitions the genome into bins, extracts the raw readcounts, filters using mappability and GC content, and normalizes the read counts using normal cells identified in the sample. It is executed as follows:
```
seacon prep_readcount -o out_dir -i bam_path -r ref.fa -b 5000000 --map-file hg38_mappability.bigWig [options]
```
If a matched-normal sample is available, use the `--normal-bam` parameter. The readcounts can also be normalized without using the normal cells with the `--no-normal` flag.

### prep_baf mode
The *prep_baf* mode initializes input precomputed BAFs, or executes the BAF estimation procedure of CHISEL using a vcf file with phased SNPs. It is executed as follows:
```
seacon prep_baf -o out_dir -i bam_path --vcf phased_snps.vcf.gz --block-size 100 [options] 
```
Note that the block size is in kb. If providing preocmputed BAFs, use the `--precomputed-baf` parameter followed by the path to the input file. The file should be in tsv format with headers chrom, start, end, cell, BAF. 

### call mode
Once the previous two commands are run, the *call* mode infers the copy numbers.
```
seacon call -o out_dir [options]
```
 
### pseudonormal mode
If no matched-normal sample is available, SEACON can construct a pseudonormal sample which can be used to construct an input vcf file. This mode should be called after running the *prep_readcount* mode, and is executed as follows:
```
seacon pseudonormal -i bam_dir -o out_dir
```
After this file is generated, use the provided `gen_vcf.sh` script in the *scripts* directory to identify the phased SNPs (see the section [Obtaining phased SNPs](#obtaining-phased-snps)). When this is complete, continue by running SEACON with the *prep_baf* mode.

### output files
After running the full SEACON pipeline, the final allele-specific copy number profiles are saved in the *calls.tsv* file, which has the following format:

cell &nbsp; chrom &nbsp; start &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; end &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; CN <br>
cell1 &nbsp; chr1	&nbsp;&nbsp;&nbsp; 1000001 &nbsp; 2000000 &nbsp;&nbsp;&nbsp; 1,1 <br>
cell1 &nbsp; chr1	&nbsp;&nbsp;&nbsp; 2000001 &nbsp; 3000000 &nbsp;&nbsp;&nbsp; 1,2 <br>
cell1 &nbsp; chr1	&nbsp;&nbsp;&nbsp; 3000001 &nbsp; 4000000 &nbsp;&nbsp;&nbsp; 2,3

The copy number profiles are also saved more compactly in the *calls_flat.tsv* file. Here, each row contains the cell name followed by its copy numbers across the genome. The *ith* column corresponds to the genomic bin specified by the *ith* line in the *filtered_regions.bed* file.

The readcounts, RDRs, and BAFs are saved in files *readcounts.tsv*, *RDR.tsv*, and *BAF.tsv*, respectively, also in compact format. Additionally, the set of normal cells identified in the sample are saved in the *diploid.txt* file.

If the BAF estimation procedure of CHISEL is used, the output of chisel can be found in file *chisel_out/calls/calls.tsv*. Please see the CHISEL documentation for more details.

## Utilities

### Raw data preprocessing
If raw reads are acquired in fastq format, we recommend the following steps (note that the bioinformatics tools *bwa* and *samtools* are required):
```
bwa mem -t 8 ref.fa cell.fastq.gz > cell.sam
samtools view -q 40 -hbS -o cell.bam cell.sam
samtools sort -o cell.sorted.bam cell.bam
samtools index cell.sorted.bam
```

### Obtaining phased SNPs
If a matched-normal sample is available, or if a pseudonormal sample is constructed using the *pseudonormal* mode of SEACON, it is necessary to obtain a set of phased SNPs to compute the BAFs. A pipeline for obtaining a vcf file with phased SNPs is provided in the script `gen_vcf.sh` located in the *scripts* directory. This script requires *bcftools* and additionally a phasing tool such as *Eagle2* or *shapeit*.

In the provided script `gen_vcf.sh`, the tool Eagle2 is used, and a prebuilt binary for linux can be installed [here](https://alkesgroup.broadinstitute.org/Eagle/). Using reference-based phasing locally requires a reference panel to be downloaded, for example [this reference](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/) from 1000 Genomes. Reference-free phasing is also possible. Additionally, Eagle2 can be used online through the the [Michigan Imputation Server](https://imputationserver.sph.umich.edu/index.html#!).

### Example commands
The following is an example workflow of using SEACON
```
# Normalize readcounts and identify normal cells
seacon prep_readcount -o seacon_output -i dataset1_bams -r ../reference/hg38.fa -b 1000000 -map-file ../reference/mappability/hg38_mappability.bigWig -P 12

# Build pseudobulk sample from the normal cells identified with the previous step. Output is located in seacon_output/pseudonormal.bam
seacon pseudonormal -o seacon_output -i dataset1_bams

# Customize the provided script to obtain phased snps from the pseudonormal.bam file
bash scripts/gen_vcf.sh

# Compute the BAFs using the phased SNPs
seacon prep_baf -o seacon_output -i dataset1_bams --vcf phased_snps.vcf.gz --block-size 100

# With the input data ready, infer the copy numbers. Here, we are using a component merge tolerance threshold of 0.2 and a breakpoint filter of 5 cells to overcome noisy measurements. 
seacon call -o seacon_output --tolerance 0.2 --upper-filter 5 --max-wgd 2
```

## Simulations
Datasets were generated using the [CNAsim](https://github.com/samsonweiner/CNAsim) simulator. Information for recreating datasets used in the paper can be found in the README.md file in the *simulation* directory.
