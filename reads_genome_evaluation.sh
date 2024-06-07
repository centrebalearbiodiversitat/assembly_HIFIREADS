# Reads quality control using NanoPlot and LongQC 
## Genome Size Estimation using KMC, kmc_tools, Genomescope2 and Smudgeplot

HIFI_READS=$1
THREADS=$2
if [ -z "$1" ]; then
    echo "Usage: $0 <path_to_hifi_reads>"
    exit 1
fi

# Quality Control
mkdir -p reads_quality
NanoPlot -t ${THREADS} --fastq ${HIFI_READS} -o reads_quality/NanoPlot_hifi
python3 /opt/LongQC/longQC.py sampleqc --ncpu ${THREADS} -o reads_quality/LongQC_hifi -x pb-hifi ${HIFI_READS} 
echo "longqc error"
cd ..

# Genome size assessment using KMC, KMC_tools, GenomeScope2 and SmudgePlot
mkdir -p genome_metrics
cd genome_metrics
mkdir -p  kmc_temp #kmc temporary directory
# kmc using 21-mers and 24 cores
kmc -k21 -m${THREADS} ${HIFI_READS} kmc_result ./kmc_temp > kmc_output

mv kmc_result.* kmc_temp
# To be able to apply GenomeScope2: kmc_tools that will produce a histogram of k-mers occurrences.
kmc_tools transform kmc_temp/kmc_result histogram reads.histo -cx10000

# GenomeScope is in mamba, initiate the mamba env from shell
source /opt/mamba/mambaforge/etc/profile.d/conda.sh

conda activate genome_scope

# input must be the reads.histo output from kmc_tools, -k is kmer length, -p is ploidy
genomescope2 -i reads.histo -k 21 -p 2 -o genomescope_out

# Use smudgeplot to check ploidy
mkdir -p smudgeplot
cp reads.histo smudgeplot/
cp kmc_temp/kmc_result.* smudgeplot/

cd smudgeplot
#1 Activate environment where smudgeplot is installed (Genome_scope in our case)

#2 Set the different L and U values to extract kmers in the coverage range from L to U using kmc_tools
L=$(smudgeplot.py cutoff reads.histo L)
U=$(smudgeplot.py cutoff reads.histo U)

kmc_tools transform kmc_result -ci"$L" -cx"$U" dump -s kmc_result_L"$L"_U"$U".dump
# run smudgeoplot hetkmers on the reduced file to compute the set of kmer pairs
smudgeplot.py hetkmers -o kmc_result_L"$L"_U"$U" < kmc_result_L"$L"_U"$U".dump
smudgeplot.py plot kmc_result_L"$L"_U"$U"_coverages.tsv
mamba deactivate


echo "Genome size estimation, ploidy estimation and analysis results are in the genome_metrics directory, reads quality stats are placed in reads_quality folder"."


