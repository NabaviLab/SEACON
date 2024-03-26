#!/bin/bash

indir=bam_dir
outdir=chisel_out
phasefile=phases.tsv
numproc=20
binsize=1Mb
blocksize=100kb
ref=hg38.fa

if ! [ -d $outdir ]
then
    mkdir $outdir
fi

chisel_prep -r $ref -o barcodes.bam -j $numproc -x $outdir $indir/*.bam
chisel_pseudonormal -r $ref -j $numproc -e 0.9 -n pseudonormal.bam -x $outdir $outdir/barcodes.bam
chisel -x $outdir -t $outdir/barcodes.bam -n $outdir/pseudonormal.bam -r $ref -l $phasefile -b $binsize -j $numproc -k $blocksize