##2 ASSEMBLY STEP AND ASSEMBLY QUALITY CHECK
## usage ./assembly_quality_check.sh "HIFI_reads" "HIC_1_reads" "HIC_2_reads" "cpu" "assembly_name"

HIFI_READS=$1
HIC_1=$2
HIC_2=$3
THREADS=$4
ASM_NAME=$5

mkdir hifiasm_HiC
cd hifiasm_HiC
hifiasm -o ${ASM_NAME}_assembly_hic_prim -t ${THREADS} --h1 ../${HIC_1} --h2 ../${HIC_2} ../${HIFI_READS} --primary

# get fasta files from hifiasm assemblies
awk '/^S/{print ">"$2;print $3}' ${ASM_NAME}_assembly_hic_prim.hic.p_ctg.gfa > ${ASM_NAME}_asm_prim.hic.p_ctg.fasta
awk '/^S/{print ">"$2;print $3}' ${ASM_NAME}_assembly_hic_prim.hic.a_ctg.gfa > ${ASM_NAME}_asm_prim.hic.a_ctg.fasta
cd ..

## Metrics using QUAST + Busco + gfastats

# 3.1 BUSCO of contigs obtained from assembler
mkdir -p assembly_metrics/busco_out
cd assembly_metrics

source /opt/mamba/mambaforge/etc/profile.d/conda.sh
conda activate busco
busco -i ../asm_prim.hic.p_ctg.fasta -o busco_out --lineage arthropoda_odb10 -c ${THR} -m geno

conda deactivate
# 3.2 QUAST metrics of the purge dups 3rd step.
mkdir -p quast
quast.py ../asm_prim.hic.p_ctg.fasta --large --est-ref-size 1180556766 -o quast

/opt/gfastats/build/bin/gfastats asm_prim.hic.p_ctg.fasta > asm_prim.hic.p_ctg.gfastats

cd ..

## Purge_dups. In order to perform purge_dups there are many steps that have to be done before.

mkdir purgedups
cd purgedups
prim_asm='../hifiasm_HiC/asm_prim.hic.p_ctg.fasta'
sec_asm='../hifiasm_HiC/asm_prim.hic.a_ctg.fasta'
# step1. split an assembly and do a self-self alignment.
split_fa ${prim_asm} > ${prim_asm}.split
minimap2 -xmap-hifi -DP ${prim_asm}.split ${prim_asm}.split | gzip -c > ${prim_asm}.split.self.paf.gz

# step2. purge haplotigs and overlaps
purge_dups -2 -T cutoffs -c PB.base.cov ${prim_asm}.split.self.paf.gz > dups.bed 2>  purge_dups.log

# step3. get purged primary and haplotig sequences from draft assembly
get_seqs -e dups.bed ${prim_asm}

mkdir -p quality_stats/busco
mkdir -p quality_stats/quast

conda activate busco
# Applying BUSCO to the first 3 steps output.
busco -i purged.fa -o quality_stats/busco --lineage arthropoda_odb10 -c ${THR} -m geno

conda deactivate

quast.py hap.fa --large --est-ref-size 1180556766 -o quality_stats/quast


