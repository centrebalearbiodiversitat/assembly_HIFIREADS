#!/bin/bash
echo "Step 5 -- Preparing Blobtoolkit"
mkdir -p BlobTool
cd BlobTool
YAHS_OUTPUT=$1
HIC_1=$2
HIC_2=$3
ASM_NAME=$4
tax_id=$5

cp ../YAHS/${YAHS_OUTPUT} .
# First,   we need a  alignment.bam, a blast_out and the taxdump database!
# Alignment:
#build index of the genome
bwa index ${YAHS_OUTPUT}
bwa mem -5SPM -T30 -t24 ${YAHS_OUTPUT} ../${HIC_1} ../${HIC_2} | samtools view -Shb -@ 6 > alignment.bam
samtools sort -@ 12 -o alignment_sorted.bam alignment.bam
samtools flagstat alignment_sorted.bam > alignment_sorted.bam.flagstats
samtools index alignment_sorted.bam
# blastn:
blastn -db nt -query ${YAHS_OUTPUT} -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 1 -max_hsps 1 -evalue 1e-25 -num_threads 64 -out ${YAHS_OUTPUT}_blast.out
# Tambien  se necesita descargar la base de datos taxdump:
#wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
#tar -xzvf taxdump.tar.gz
# RUN BLOBTOOLKIT:
# Step 1. Create blobtools database:
blobtools create --fasta ${YAHS_OUTPUT}  --meta meta.yaml --taxid $tax_id --taxdump . ${ASM_NAME}
# Step 2. Add blast hits
blobtools add  --hits ${YAHS_OUTPUT}_blast.out --taxrule bestsumorder --taxdump . ${ASM_NAME}
# Step 3. Add coverage
# samtools index -c alignment_sorted.bam ## looks for a .bam.csi
blobtools add --cov alignment_sorted.bam ${ASM_NAME}
# Step 4. Add busco scores
blobtools add --busco ../YAHS/busco_out/run_*/full_table.tsv ${ASM_NAME}
# Step 5. Open dataset in Blobtoolkit viewer
blobtools host ${ASM_NAME}

echo "Step 5 -- DONE. check http:localhost:8080"
