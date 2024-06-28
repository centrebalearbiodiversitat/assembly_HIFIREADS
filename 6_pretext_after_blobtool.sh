#!/bin/bash

FASTA=$1 #Blobtool output
HIC_1=$2
HIC_2=$3
#Create an index from the YAHS output
bowtie2-build ${FASTA} -o ${FASTA}_bowtie
# Align reads to reference genome
bowtie2 -x ${FASTA}_bowtie -1 ../${HIC_1} -2 ../${HIC_2} -S ${FASTA}_aligned.sam
# Convert SAM to BAM and sort
samtools view -bS ${FASTA}_aligned.sam | samtools sort -o ${FASTA}_alignment_sorted.bam -
samtools flagstat ${FASTA}_alignment_sorted.bam > ${FASTA}_align_sorted.bam.flagstats

# PretextMap needs an alignment .bam.
# Obtain .pretext:
samtools view -h ${FASTA}_aligment.sorted.bam | /opt/PretextMap/builddir/PretextMap -o map.pretext --sortby length --sortorder descend --mapq 10
#Pretext coverage 
# Index BAM file
samtools index ${FASTA}_alignment_sorted.bam
# Compute coverage
bedtools genomecov -ibam ${FASTA}_alignment_sorted.bam -bg > ${FASTA}_coverage.bedgraph
# Add coverage.graph to the Map
awk '$5=="U" {print $1 "\t" $2 "\t" $3 "\t1000" }' ${BLOBTOOL_AGP}.agp > gaps.bed
cat gaps.bed | /opt/PretextGraph/bin/PretextGraph addGaps -i map.pretext -n "Gaps"
cat coverage.bedgraph | PretextGraph -i map.pretext -n "coverageMap"
cat telomeres.bed | PretextGraph -i map.pretext  -n "Telomeres"
#Open PretextView and the pretext.map

