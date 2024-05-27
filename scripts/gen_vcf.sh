#!/bin/bash

outdir=out1
bampath=out1/pseudonormal.bam
reference=../ref/hg38.fa
snpdir=$outdir/snp
panel=../reference/1000G/panel

if [ ! -d $snpdir ]
then
    mkdir $snpdir
fi

# Call variants
bcftools mpileup --skip-indels --ignore-RG -f $reference $bampath | bcftools call -mv -Ob -o $snpdir/calls.bcf
bcftools index $snpdir/calls.bcf

for chr in {1..22}; do
    bcftools view -Ob -r chr${chr} -o snps/chr${chr}.bcf $snpdir/calls.bcf
    bcftools index $snpdir/chr${chr}.bcf
done

# Phase with eagle2
for chr in {1..22}; do
    ./eagle --vcfTarget $snpdir/snps/chr${chr}.bcf --vcfRef $panel/chr${chr}.bcf.gz --geneticMapFile ../tools/Eagle2/tables/genetic_map_hg38_withX.txt.gz --outPrefix $snpdir/chr${chr}_phased --vcfOutFormat b
done

#Concatenate all chromosomes into a single file.
bcftools concat -Ob -o all_phased.bcf phased_snps/chr{1..22}.bcf

#Filter for heterozygous sites
bcftools view -Ob -g het -o all_phased_het.bcf all_phased.bcf
bcftools view -Oz -o all_phased_het.vcf.gz all_phased_het.bcf