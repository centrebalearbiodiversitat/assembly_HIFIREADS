#! /bin/bash
## Pipeline created by CBB to assemble WHOLE GENOME with HiFi PACBIO long-reads and HiC.
## it would be necessary to change the lineage for BUSCO and the genome size after using genome_scope.
# example full_pipeline.sh /reads/hifi.fastq hic_1.fq hic_2.fq 30 species
#first input must be HiFi reads
#second input. HIC_1.fq
#third input. HIC_2.fq
#fourth input must be threads used
#third input assembly must be "assembly name being for example "Raja_radula" "

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
MAIN_DIR=$(pwd)
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

genome_size=$(awk '/Genome Haploid Length/ { gsub(",", "", $4); gsub(",", "", $6); print ($4 + $6) / 2 }' genomescope_out/summary.txt)

# Use smudgeplot to check ploidy
mkdir -p smudgeplot
cd smudgeplot
export PATH=$PATH:/opt/FASTK/:/opt/smudgeplot-sploidyplot/exec/
FastK -T4 -k21 -t4 -M16 ../../../${HIFI_READS} -NFastK_table
PloidyPlot -o${ASM_NAME} FastK_table
smudgeplot.py plot -t ${ASM_NAME} -o smudge ${ASM_NAME}_text.smu

mamba deactivate


echo "Step 1 - DONE"
echo "Genome size estimation, ploidy estimation and analysis results are in the genome_metrics directory, reads quality stats are placed in reads_quality folder"

echo "Starting step 2 - assembly"
##2 ASSEMBLY STEP AND ASSEMBLY QUALITY CHECK

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
busco -i ../hifiasm_HiC/${ASM_NAME}_asm_prim.hic.p_ctg.fasta -o busco_out --auto-lineage-euk -c ${THREADS} -m geno -f

## obtain the best lineage for further steps
lineage=$(ls busco_out | grep 'short_summary.specific.' | grep '.txt'| cut -d'.' -f3,3)

conda deactivate
# 3.2 QUAST metrics of the purge dups 3rd step.
mkdir -p quast
quast.py ../hifiasm_HiC/${ASM_NAME}_asm_prim.hic.p_ctg.fasta --large --est-ref-size ${genome_size} -o quast

/opt/gfastats/build/bin/gfastats ../hifiasm_HiC/${ASM_NAME}_asm_prim.hic.p_ctg.fasta > ${ASM_NAME}_asm_prim.hic.p_ctg.gfastats

cd ..

## Purge_dups. In order to perform purge_dups there are many steps that have to be done before.

mkdir purgedups
cd purgedups
cp ../hifiasm_HiC/${ASM_NAME}_asm_prim.hic.p_ctg.fasta .
prim_asm='../hifiasm_HiC/${ASM_NAME}_asm_prim.hic.p_ctg.fasta'
# step1. split an assembly and do a self-self alignment.
split_fa ${prim_asm} > ${prim_asm}.split
minimap2 -xmap-hifi -DP ${prim_asm}.split ${prim_asm}.split | gzip -c > ${prim_asm}.split.self.paf.gz

# step2. purge haplotigs and overlaps
purge_dups -2 -T cutoffs -c PB.base.cov ${prim_asm}.split.self.paf.gz > dups.bed 2>  purge_dups.log

# step3. get purged primary and haplotig sequences from draft assembly
get_seqs -e dups.bed ${prim_asm}

mkdir -p quality_stats/busco_out
mkdir -p quality_stats/quast

conda activate busco
# Applying BUSCO to the first 3 steps output.
busco -i purged.fa -o quality_stats/busco_out --lineage ${lineage} -c ${THREADS} -m geno -f

conda deactivate

quast.py purged.fa --large --est-ref-size ${genome_size} -o quality_stats/quast

/opt/gfastats/build/bin/gfastats purged.fa > ${ASM_NAME}_purged.gfastats

echo "Step 2 -- Done. Assembly, purge duplicates and quality assessment prior and after the following steps have been performed."

echo "Step 3 -- Decontamination before scaffolding"

cd ${MAIN_DIR}
mkdir -p decontamination/fcs_output
mkdir -p decontamination/whokaryote_output

## obtain tax_id 
tax_id=$(python3 get_taxon_id.py "${ASM_NAME}")

python3 /opt/fcs/fcs.py screen genome --fasta purgedups/purged.fa --out-dir decontamination/fcs_output/ --gx-db /opt/fcs/gxdb/ --tax-id ${tax_id}
# Delete contaminants:
cat purgedups/purged.fa | python3 /opt/fcs/fcs.py clean genome --action-report decontamination/fcs_output/purged.${tax_id}.fcs_gx_report.txt --output /decontamination/fcs_output/${ASM_NAME}_FCS_clean.fasta --contam-fasta-out decontamination/fcs_output/${ASM_NAME}_FCS_contam.fasta

source /opt/mamba/mambaforge/etc/profile.d/conda.sh
mamba activate whykaryote
# Whokaryote
whokaryote.py --contigs decontamination/whokaryote_output/${ASM_NAME}_FCS_clean.fasta --outdir decontamination/whokaryote_output --f --minsize 10000 --model T
mamba deactivate
cd ${MAIN_DIR}
echo "Step 3 -- Done. Decontamination prior scaffolding has been performed"

echo "Step 4 -- Scaffolding using YAHS"

mkdir bwa_output
cd bwa_output

cp ../decontamination/whokaryote_output/eukaryotes.fa .
file=eukaryotes.fa
bwa index ${file}
samtools faidx $file
## checking quality with fastqc for HI-C reads
fastqc ../${HIC_1} ../${HIC_2}
# use BWA-MEM to align the Hi-C paired-end reads to reference sequences
# Step1. HiC reads from FASTQC to BAM
HIC1_basename=$(basename "${HIC_1%.*}")
HIC2_basename=$(basename "${HIC_2%.*}")
HIC1_output="${HIC1_basename}.filtered.bam"
HIC2_output="${HIC2_basename}.filtered.bam"

bwa mem -t ${THREADS} ${file} ../${HIC_1} | samtools view -@ ${THREADS} -Sb - > "${HIC1_basename}.bam"
bwa mem -t ${THREADS} ${file} ../${HIC_2} | samtools view -@ ${THREADS} -Sb - > "${HIC2_basename}.bam"

#Step2. Retain only the portion of the chimeric read that maps in the 5'-orientation in relation to its read orientation.

samtools view -h "${HIC1_basename}.bam" | perl /opt/scripts/filter_five_end.pl | samtools view -Sb > "${HIC1_output}"
samtools view -h "${HIC2_basename}.bam" | perl /opt/scripts/filter_five_end.pl | samtools view -Sb > "${HIC2_output}"

# Now we pair the filtered single-end Hi-C using "two_read_bam_combiner.pl"

REF='${file}'
FAIDX='$REF.fai'
perl /opt/scripts/two_read_bam_combiner.pl "${HIC1_output}" "${HIC2_output}"  samtools 10 | samtools view -bS -t $FAIDX | samtools sort -@ ${THREADS} -o HiC1_HiC2_combined.bam

YAHS_OUTPUT=$(yahs eukaryotes.hic.asm HiC1_HiC2_combined.bam)

source /opt/mamba/mambaforge/etc/profile.d/conda.sh
conda activate busco
busco -i  ${YAHS_OUTPUT} -c ${THREADS} --lineage ${lineage} -o busco_out -m geno 
conda deactivate
/opt/gfastats/build/bin/gfastats -i ${YAHS_OUTPUT} > ${YAHS_OUTPUT}.gfastats

echo "Step 4 -- DONE"

echo "Step 5 -- Preparing Blobtoolkit"
                                                                                                         
# First,   we need a  alignment.bam, a blast_out and the taxdump database!
# Alignment:
#build index of the genome
bwa index $YAHS_OUTPUT
bwa mem -5SPM -T30 -t24 $YAHS_OUTPUT ${HIC_1} ${HIC_2} | samtools view -Shb -@ 6 > alignment.bam
samtools sort -@ 12 -o alignment_sorted.bam alignment.bam
samtools flagstat alignment_sorted.bam > alignment_sorted.bam.flagstats
samtools index alignment_sorted.bam
# blastn:
blastn -db nt -query $outfile_fasta -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 1 -max_hsps 1 -evalue 1e-25 -num_threads 64 -out blast.out
# Tambien  se necesita descargar la base de datos taxdump:
#wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
#tar -xzvf taxdump.tar.gz
# RUN BLOBTOOLKIT:
# Step 1. Create blobtools database:
blobtools create --fasta $YAHS_OUTPUT  --meta meta.yaml --taxid $tax_id --taxdump . ${ASM_NAME}
# Step 2. Add blast hits
blobtools add  --hits blast.out --taxrule bestsumorder --taxdump . ${ASM_NAME}
# Step 3. Add coverage
# samtools index -c alignment_sorted.bam ## looks for a .bam.csi
blobtools add --cov alignment_sorted.bam ${ASM_NAME}
# Step 4. Add busco scores
blobtools add --busco busco_out/run_*/full_table.tsv ${ASM_NAME}
# Step 5. Open dataset in Blobtoolkit viewer
blobtools host `pwd`


