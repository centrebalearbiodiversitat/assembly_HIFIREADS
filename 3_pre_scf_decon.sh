#!/bin/bash
# Use ./decontamination.sh SPECIES_NAME
echo "Step 3 -- Decontamination before scaffolding"
ASM_NAME=$1
MAIN_DIR=${ASM_NAME}
cd ${MAIN_DIR}
mkdir -p decontamination/fcs_output
mkdir -p decontamination/whokaryote_output

## obtain tax_id 
tax_id=$(python3 get_taxon_id.py "${ASM_NAME}")

python3 /opt/fcs/dist/fcs.py screen genome --fasta purgedups/purged.fa --out-dir decontamination/fcs_output/ --gx-db /opt/fcs/gxdb/ --tax-id ${tax_id}

# Delete contaminants:
cat purgedups/purged.fa | python3 /opt/fcs/fcs.py clean genome --action-report decontamination/fcs_output/fcs_output/purged.${tax_id}.fcs_gx_report.txt --output decontamination/fcs_output/${ASM_NAME}_FCS_clean.fasta --contam-fasta-out decontamination/fcs_output/${ASM_NAME}_FCS_contam.fasta

source /opt/mamba/mambaforge/etc/profile.d/conda.sh
conda activate whokaryote
# Whokaryote
whokaryote.py --contigs decontamination/fcs_output/${ASM_NAME}_FCS_clean.fasta --outdir decontamination/whokaryote_output --f --minsize 10000 --model T
conda deactivate
cd ${MAIN_DIR}
echo "Step 3 -- Done. Decontamination prior scaffolding has been performed"



