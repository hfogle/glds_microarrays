#!/bin/bash


workdir=pwd
tmpdir=$(mktemp -d -t ci-XXXXXXXXXXXXXXXXXXXX)
cd $tmpdir
echo $tmpdir
#eval "$(conda shell.bash hook)"
#conda activate glds_microarrays

retrieve_isa_from_genelab.py --accession $1 --to-Microarray-runsheet --allow-missing-columns 

echo $1
cd $workdir
Rscript --vanilla glds_microarrays.R --staging $tmpdir/*.csv --glds $1 --isa $tmpdir/*ISA.zip

rm -rf $tmpdir