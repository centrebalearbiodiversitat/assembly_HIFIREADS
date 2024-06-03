! /bin/bash

## Genome Size Estimation using KMC, kmc_tools, Genomescope2 and Smudgeplot
## place in a folder names "reads". 
READS = $1

if [ -z "$1" ]; then
    echo "Usage: $0 <path_to_reads>"
    exit 1
fi

mkdir -p genome_metrics
cd genome_metrics
##step 1. kmer and genome size estimation using KMC, kmc_tools, genomescope and smudgeplot
mkdir -p  kmc_temp #kmc temporary directory
# kmc using 21-mers and 24 cores
kmc -k21 -m24 $READS kmc_result ./kmc_temp > kmc_output

mv kmc_result.* kmc_temp
# To be able to apply GenomeScope2: kmc_tools that will produce a histogram of k-mers occurrences.
kmc_tools transform kmc_temp/kmc_result histogram reads.histo -cx10000

# GenomeScope is in mamba, initiate the mamba env from shell

source /opt/mamba/mambaforge/etc/profile.d/conda.sh
mamba activate genome_scope

# input must be the reads.histo output from kmc_tools, -k is kmer length, -p is ploidy
genomescope2 -i reads.histo -k 21 -p 2 -o genome_out

# Use smudgeplot to check ploidy
mkdir -p smudgeplot
cp genomescope_out/reads.histo smudgeplot/
cp kmc_output/kmc_out.* smudgeplot/

cd smudgeplot
#1 Activate environment where smudgeplot is installed (Genome_scope in our case)

#2 Set the different L and U values to extract kmers in the coverage range from L to U using kmc_tools
L=$(smudgeplot.py cutoff reads.histo L)
U=$(smudgeplot.py cutoff reads.histo U)

kmc_tools transform kmc_out -ci"$L" -cx"$U" reduce kmc_out_L"$L"_U"$U"
# run smudge_pairs on the reduced file to compute the set of kmer pairs
smudge_pairs kmc_out_L"$L"_U"$U" kmc_out_L"$L"_U"$U"_coverages.tsv kmc_out_L"$L"_U"$U"_pairs.tsv > kmc_out_L"$L"_U"$U"_familysizes.tsv
smudgeplot.py plot kmc_out_L"$L"_U"$U"_coverages.tsv

mamba deactivate

echo "Genome size estimation, ploidy estimation and analysis results are in the genome_metrics directory."
