#!/bin/bash

# FCS- NCBI
# Screen the genome:
python3 /opt/fcs/fcs.py screen genome --fasta purged.fa --out-dir ./output/ --gx-db /opt/fcs/gxdb/ --tax-id 203899
# Delete contaminants:
cat purged.fa | python3 /opt/fcs/fcs.py clean genome --action-report ./output/purged.203899.fcs_gx_report.txt --output clean.fasta --contam-fasta-out contam.fasta

# Whokaryote
whokaryote.py --contigs clean.fasta --outdir whokaryote_output --f --minsize 10000 --model T
