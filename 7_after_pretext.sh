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


