#! /bin/bash
## Pipeline created by CBB to assemble WHOLE GENOME with HiFi PACBIO long-reads and HiC.
## it would be necessary to change the lineage for BUSCO and the genome size after using genome_scope.
# example full_pipeline.sh /reads/hifi.fastq hic_1.fq hic_2.fq 30 species
#first input must be HiFi reads
#second input. HIC_1.fq
#third input. HIC_2.fq
#fourth input must be threads used
#third input assembly must be "assembly name"

echo "Starting pipeline"
## 1. reads_genome_evaluation
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
#mkdir -p reads_quality
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
awk '/Genome Haploid Length/ { gsub(",", "", $4); gsub(",", "", $6); print ($4 + $6) / 2 }' genomescope_out/summary.txt > genomescope_out/estimation_length.txt

genome_size=$(cat estimation_length.txt)

# Use smudgeplot to check ploidy
mkdir -p smudgeplot
cd smudgeplot
#1 Activate environment where smudgeplot is installed (Genome_scope in our case)

#2 Set the different L and U values to extract kmers in the coverage range from L to U using kmc_tools
export PATH=$PATH:/opt/FASTK/:/opt/smudgeplot-sploidyplot/exec/
FastK -T4 -k21 -t4 -M16 ../../${HIFI_READS} -NFastK_table
PloidyPlot -o${ASM_NAME} FastK_table
smudgeplot.py plot -t ${ASM_NAME} -o smudge ${ASM_NAME}_text.smu

#mamba deactivate


echo "Step 1 - DONE"
echo "Genome size estimation, ploidy estimation and analysis results are in the genome_metrics directory, reads quality stats are placed in reads_quality folder"

echo "Starting step 2 - assembly"
##2 ASSEMBLY STEP AND ASSEMBLY QUALITY CHECK
## usage ./assembly_quality_check.sh "HIFI_reads" "HIC_1_reads" "HIC_2_reads" "cpu"

mkdir hifiasm_HiC
cd hifiasm_HiC
hifiasm -o ${ASM_NAME}_assembly_hic_prim -t ${THREADS} --h1 ../${HIC_1} --h2 ../${HIC_2} ../${HIFI_READS} --primary

# get fasta files from hifiasm assemblies
awk '/^S/{print ">"$2;print $3}' ${ASM_NAME}_assembly_hic_prim.hic.p_ctg.gfa > ${ASM_NAME}_asm_prim.hic.p_ctg.fasta
awk '/^S/{print ">"$2;print $3}' ${ASM_NAME}_assembly_hic_prim.hic.a_ctg.gfa > ${ASM_NAME}_asm_prim.hic.a_ctg.fasta
cd ..

## Metrics using QUAST + Busco + gfastats

# 3.1 BUSCO of contigs obtained from assembler
mkdir -p assembly_metrics/busco_out
cd assembly_metrics

source /opt/mamba/mambaforge/etc/profile.d/conda.sh
conda activate busco
busco -i ../${ASM_NAME}_asm_prim.hic.p_ctg.fasta -o busco_out --lineage arthropoda_odb10 -c ${THREADS} -m geno

conda deactivate
# 3.2 QUAST metrics of the purge dups 3rd step.
mkdir -p quast
quast.py ../${ASM_NAME}_asm_prim.hic.p_ctg.fasta --large --est-ref-size ${genome_size} -o quast

/opt/gfastats/build/bin/gfastats ${ASM_NAME}_asm_prim.hic.p_ctg.fasta > ${ASM_NAME}_asm_prim.hic.p_ctg.gfastats

cd ..

## Purge_dups. In order to perform purge_dups there are many steps that have to be done before.

mkdir purgedups
cd purgedups
prim_asm='../hifiasm_HiC/${ASM_NAME}_asm_prim.hic.p_ctg.fasta'
sec_asm='../hifiasm_HiC/${ASM_NAME}_asm_prim.hic.a_ctg.fasta'
# step1. split an assembly and do a self-self alignment.
split_fa ${prim_asm} > ${prim_asm}.split
minimap2 -xmap-hifi -DP ${prim_asm}.split ${prim_asm}.split | gzip -c > ${prim_asm}.split.self.paf.gz

# step2. purge haplotigs and overlaps
purge_dups -2 -T cutoffs -c PB.base.cov ${prim_asm}.split.self.paf.gz > dups.bed 2>  purge_dups.log

# step3. get purged primary and haplotig sequences from draft assembly
get_seqs -e dups.bed ${prim_asm}

mkdir -p quality_stats/busco
mkdir -p quality_stats/quast

conda activate busco
# Applying BUSCO to the first 3 steps output.
busco -i purged.fa -o quality_stats/busco --lineage arthropoda_odb10 -c ${THREADS} -m geno

conda deactivate

quast.py hap.fa --large --est-ref-size ${genome_size} -o quality_stats/quast
