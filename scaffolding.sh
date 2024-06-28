#!/bin/bash
echo "Step 4 - scaffolding"
purged=$1
HIC1=$2
HIC2=$3
THREADS=$4
mkdir bwa_output
cd bwa_output

cp ../purgedups/${purged} .
bwa index ${purged}
samtools faidx ${purged}
## checking quality with fastqc for HI-C reads
fastqc ../${HIC1} ../${HIC2}
# use BWA-MEM to align the Hi-C paired-end reads to reference sequences
# Step1. HiC reads from FASTQC to BAM
HIC1_basename=$(basename "${HIC1%.*}")
HIC2_basename=$(basename "${HIC2%.*}")
HIC1_output="${HIC1_basename}.filtered.bam"
HIC2_output="${HIC2_basename}.filtered.bam"

bwa mem -t ${THREADS} ${purged} ../${HIC1} | samtools view -@ ${THREADS} -Sb - > "${HIC1_basename}.bam"
bwa mem -t ${THREADS} ${purged} ../${HIC2} | samtools view -@ ${THREADS} -Sb - > "${HIC2_basename}.bam"

#Step2. Retain only the portion of the chimeric read that maps in the 5'-orientation in relation to its read orientation.

samtools view -h "${HIC1_basename}.bam" | perl /opt/scripts/filter_five_end.pl | samtools view -Sb > "${HIC1_output}"
samtools view -h "${HIC2_basename}.bam" | perl /opt/scripts/filter_five_end.pl | samtools view -Sb > "${HIC2_output}"

# Now we pair the filtered single-end Hi-C using "two_read_bam_combiner.pl"

REF='/home/bioinfo/tethysbaena_scabra_assembly/bwa_output/purged.fa'
FAIDX='$REF.fai'
perl /opt/scripts/two_read_bam_combiner.pl "${HIC1_output}" "${HIC2_output}"  samtools 10 | samtools view -bS -t $FAIDX | samtools sort -@ ${THREADS} -o HiC1_HiC2_combined.bam

YAHS_OUTPUT = /opt/yahs/yahs purged.hic.asm HiC1_HiC2_combined.bam

## checking scaffolding stats
mkdir -p busco_yahs
mamba activate busco
busco -i $YAHS_OUTPUT -o busco_yahs -m geno --lineage ${lineage} -c ${THREADS} -f 
mamba deactivate

/opt/gfastats/build/gfastats ${YAHS_OUTPUT} > ${YAHS_OUTPUT}.gfastats

echo "Step 4 -- Scaffolding DONE"

