#!/bin/bash

FASTA=$1 #Blobtool output
HIC_1=$2
HIC_2=$3
#Create an index from the btk output
bwa index ${FASTA}
#bowtie2-build ${FASTA} -o ${FASTA}_bowtie
# Align reads to reference genome
#bowtie2 -x ${FASTA}_bowtie -1 ../${HIC_1} -2 ../${HIC_2} -S ${FASTA}_aligned.sam
bwa mem -5SPM -T30 -t24 ${FASTA} ../${HIC_1} ../${HIC_2} | samtools view -Shb -@ 6 > ${FASTA}_alignment.bam
samtools sort -@ 12 -o ${FASTA}_alignment_sorted.bam alignment.bam

# Convert SAM to BAM and sort
#samtools view -bS ${FASTA}_aligned.sam | samtools sort -o ${FASTA}_alignment_sorted.bam -
samtools flagstat ${FASTA}_alignment_sorted.bam > ${FASTA}_align_sorted.bam.flagstats

bwa index assembly.fasta
bwa mem -5SPM -T30 -t24 assembly.fasta ../../../../reads/Tscabra_allHiC_1.fastq.gz ../../../../reads/Tscabra_allHiC_2.fastq.gz | samtools view -Shb -@ 6 > alignment.bam
samtools sort -@ 12 -o alignment_sorted.bam alignment.bam
samtools index alignment_sorted.bam

# PretextMap needs an alignment .bam.
# Obtain .pretext:
samtools view -h ${FASTA}_aligment.sorted.bam | /opt/PretextMap/builddir/PretextMap -o map.pretext --sortby length --sortorder descend --mapq 10
#Pretext coverage 
# Index BAM file
samtools index ${FASTA}_alignment_sorted.bam
# Compute coverage
# Add coverage.graph to the Map

bedtools genomecov -ibam ${FASTA}_alignment_sorted.bam -bg > ${FASTA}_coverage.bedgraph
cat coverage.bedgraph | PretextGraph -i map.pretext -n "coverageMap"
# Add gaps.bed to the map
perl FastaToAPGFileFormatConverter.pl -i ${FASTA} -o ${FASTA}.agp
awk '$5=="U" {print $1 "\t" $2 "\t" $3 "\t1000" }' ${FASTA}.agp > gaps.bed
cat gaps.bed | /opt/PretextGraph/bin/PretextGraph addGaps -i map.pretext -n "Gaps"
#cat telomeres.bed | PretextGraph -i map.pretext  -n "Telomeres"
#Open PretextView and the pretext.map

