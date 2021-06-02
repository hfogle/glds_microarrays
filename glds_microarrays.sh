#!/bin/bash


workdir=$(pwd)
tmpdir=$(mktemp -d -t ci-XXXXXXXXXXXXXXXXXXXX)
mkdir -p $workdir/Processed_Data
mkdir -p $workdir/Processed_Data/$1
outdir=$workdir/Processed_Data/$1

echo $workdir
echo $outdir
#eval "$(conda shell.bash hook)"
#conda activate glds_microarrays
cd $outdir
retrieve_isa_from_genelab.py --accession $1 --to-Microarray-runsheet --allow-missing-columns 

echo $1
cd $workdir
Rscript --vanilla $workdir/glds_microarrays.R --staging $outdir/*.csv --glds $1 --isa $outdir/*ISA.zip --reports > $1.log 2>&1

rm -rf $tmpdir