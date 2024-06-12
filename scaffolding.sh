#!/bin/bash

mkdir bwa_output
cd bwa_output

bwa index purged.fa
## checking quality with fastqc for HI-C reads
fastqc ../reads/Tscabra_allHiC_*
# use BWA-MEM to align the Hi-C paired-end reads to reference sequences
# Step1. HiC reads from FASTQC to BAM
bwa mem -t 36 purged.fa ../reads/Tscabra_allHiC_1.fastq.gz | samtools view -@ 36 -Sb - > ../reads/Tscabra_allHiC_1.bam

bwa mem -t 36 purged.fa ../reads/Tscabra_allHiC_2.fastq.gz | samtools view -@ 36 -Sb - > ../reads/Tscabra_allHiC_2.bam

#Step2. Retain only the portion of the chimeric read that maps in the 5'-orientation in relation to its read orientation.

samtools view -h ../reads/Tscabra_allHiC_1.bam | perl /opt/filter_five_end.pl | samtools view -Sb > Tscabra_allHiC_1.filtered.bam
samtools view -h ../reads/Tscabra_allHiC_2.bam | perl /opt/filter_five_end.pl | samtools view -Sb > Tscabra_allHiC_2.filtered.bam

# Now we pair the filtered single-end Hi-C using "two_read_bam_combiner.pl"
mkdir temp_dir
REF='/home/bioinfo/tethysbaena_scabra_assembly/bwa_output/purged.fa'
FAIDX='$REF.fai'
perl /opt/two_read_bam_combiner.pl Tscabra_allHiC_1.filtered.bam Tscabra_allHiC_2.filtered.bam samtools 10 | samtools view -bS -t $FAIDX | samtools sort -@ 32 -o temp_dir/>

yahs purged.hic.asm Tscabra_HiC.bam

