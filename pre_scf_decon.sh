#!/bin/bash
# Use ./decontamination SPECIES_NAME
echo "Step 3 -- Decontamination before scaffolding"
ASM_NAME=$1
MAIN_DIR=${ASM_NAME}
cd ${MAIN_DIR}
mkdir -p decontamination/fcs_output
mkdir -p decontamination/whokaryote_output

## obtain tax_id 
tax_id=$(python3 get_taxon_id.py "${ASM_NAME}")

python3 /opt/fcs/fcs.py screen genome --fasta purgedups/purged.fa --out-dir ./fcs_output/ --gx-db /opt/fcs/gxdb/ --tax-id ${tax_id}
# Delete contaminants:
cat purged.fa | python3 /opt/fcs/fcs.py clean genome --action-report ./fcs_output/purged.${tax_id}.fcs_gx_report.txt --output ${ASM_NAME}_FCS_clean.fasta --contam-fasta-out ${ASM_NAME}_FCS_contam.fasta

# Whokaryote
whokaryote.py --contigs ${ASM_NAME}_FCS_clean.fasta --outdir whokaryote_output --f --minsize 10000 --model T

cd ${MAIN_DIR}
echo "Step 3 -- Done. Decontamination prior scaffolding has been performed"
