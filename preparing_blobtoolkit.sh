#!/bin/bash                                                                                                                                   
# From PRETEXT output 
FINAL_ASSEMBLY_TPF=$1
FINAL_ASSEMBLY_FASTA=$2
PRETEXT=$3
HIC1=$4
HIC2=$5
outfile_fasta='joined.fasta'

source agp-tpf-utils/venv/bin/activate

python agp-tpf-utils/src/tola/assembly/scripts/pretext_to_tpf.py -a ${FINAL_ASSEMBLY_TPF} -p ${PRETEXT} -o blobtool.tpf

perl /opt/rapid_join.pl -fa ${FINAL_ASSEMBLY_FASTA} -tpf blobtool.tpf -csv chrom.csv -o $outfile_fasta

# First,   we need a  alignment.bam, a blast_out and the taxdump database!
# Alignment:
#build index of the genome
bwa index $outfile_fasta
bwa mem -5SPM -T30 -t24 $outfile_fasta ${HIC1} ${HIC2} | samtools view -Shb -@ 6 > alignment.bam
samtools sort -@ 12 -o alignment_sorted.bam alignment.bam
samtools index alignment_sorted.bam
# blastn:
blastn -db nt -query $outfile_fasta -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 1 -max_hsps 1 -evalue 1e-25 -num_threads 64 -out blast.out
# Tambien  se necesita descargar la base de datos taxdump:
#wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
#tar -xzvf taxdump.tar.gz
# RUN BLOBTOOLKIT:
# Step 1. Create blobtools database:
blobtools create --fasta Bin-8.fasta  --meta meta.yaml --taxid 203899 --taxdump . Tethysbaena_scabra_assembly
# Step 2. Add blast hits
blobtools add  --hits blast.out --taxrule bestsumorder --taxdump . Tethysbaena_scabra_assembly
# Step 3. Add coverage
blobtools add --cov alignment_sorted.bam Tethysbaena_scabra_assembly
# Step 4. Add busco scores
blobtools add --busco full_table.tsv Tethysbaena_scabra_assembly
# Step 5. Open dataset in Blobtoolkit viewer
blobtools host `pwd`


