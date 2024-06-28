#! /bin/bash
## Pipeline created by CBB to assemble WHOLE GENOME with HiFi PACBIO long-reads and HiC.
## it would be necessary to change the lineage for BUSCO and the genome size after using genome_scope.
# example full_pipeline.sh /reads/hifi.fastq hic_1.fq hic_2.fq 30 species
#first input must be HiFi reads
#second input. HIC_1.fq
#third input. HIC_2.fq
#fourth input must be threads used
#fifth input assembly must be "species_name""

echo "Starting pipeline"
## Step 1. reads assess. Reads and genome assessment
echo "Starting step 1. Reads quality control, genome size and ploidy estimation"
## Reads quality control using NanoPlot and LongQC 
## Genome Size Estimation using KMC, kmc_tools, FastK, Genomescope2 and Smudgeplot

HIFI_READS=$1
HIC_1=$2
HIC_2=$3
THREADS=$4
ASM_NAME=$5
if [ -z "$1" ]; then
    echo "Usage: $0 <path_to_hifi_reads>"
    exit 1
fi

mkdir -p ${ASM_NAME}
cd ${ASM_NAME}
# Quality Control
mkdir -p reads_quality
NanoPlot -t ${THREADS} --fastq ../${HIFI_READS} -o reads_quality/NanoPlot_hifi
python3 /opt/LongQC/longQC.py sampleqc --ncpu ${THREADS} -o reads_quality/LongQC_hifi -x pb-hifi ../${HIFI_READS}

# Genome size assessment using KMC, KMC_tools, GenomeScope2 and SmudgePlot
mkdir -p genome_metrics
cd genome_metrics
mkdir -p  kmc_temp #kmc temporary directory
#kmc using 21-mers and 24 cores
kmc -k21 -m${THREADS} ../../${HIFI_READS} kmc_result ./kmc_temp > kmc_output

mv kmc_result.* kmc_temp
# To be able to apply GenomeScope2: kmc_tools that will produce a histogram of k-mers occurrences.
kmc_tools transform kmc_temp/kmc_result histogram reads.histo -cx10000

# GenomeScope is in mamba, initiate the mamba env from shell
source /opt/mamba/mambaforge/etc/profile.d/conda.sh

conda activate genome_scope

# input must be the reads.histo output from kmc_tools, -k is kmer length, -p is ploidy
genomescope2 -i reads.histo -k 21 -p 2 -o genomescope_out
genome_size=$(awk '/Genome Haploid Length/ { gsub(",", "", $4); gsub(",", "", $6); print ($4 + $6) / 2 }' genomescope_out/summary.txt)

# Use smudgeplot to check ploidy
mkdir -p smudgeplot
cd smudgeplot

export PATH=$PATH:/opt/FASTK/:/opt/smudgeplot-sploidyplot/exec/
FastK -T${THREADS}-k21 -t4 -M16 ../../../${HIFI_READS} -NFastK_table
PloidyPlot -o${ASM_NAME} FastK_table
smudgeplot.py plot -t ${ASM_NAME} -o smudge ${ASM_NAME}_text.smu

#mamba deactivate


echo "Step 1 - DONE"
echo "Genome size estimation, ploidy estimation and analysis results are in the genome_metrics directory, reads quality stats are placed in reads_quality folder"
