##2 ASSEMBLY STEP AND ASSEMBLY QUALITY CHECK
## usage ./assembly_quality_check.sh "HIFI_reads" "HIC_1_reads" "HIC_2_reads" "cpu" "assembly_name"
## check the number of expected haplotypes
HIFI_READS=$1
HIC_1=$2
HIC_2=$3
THREADS=$4
ASM_NAME=$5

mkdir hifiasm_HiC
cd hifiasm_HiC
hifiasm -o ${ASM_NAME}_assembly_hic_prim -t ${THREADS} --h1 ../${HIC_1} --h2 ../${HIC_2} ../${HIFI_READS} --primary --n-hap 40

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
busco -i ../hifiasm_HiC/${ASM_NAME}_asm_prim.hic.p_ctg.fasta -o busco_out --auto-lineage-euk -c ${THREADS} -m geno

## obtain the best lineage for further steps
lineage=$(ls busco_out | grep 'short_summary.specific.' | grep '.txt'| cut -d'.' -f3,3)

conda deactivate
# 3.2 QUAST metrics of the purge dups 3rd step.
mkdir -p quast
genome_size=220000000
quast.py ./hifiasm_HiC/${ASM_NAME}_asm_prim.hic.p_ctg.fasta --large --est-ref-size $genome_size -o quast

/opt/gfastats/build/bin/gfastats ./hifiasm_HiC/${ASM_NAME}_asm_prim.hic.p_ctg.fasta > asm_prim.hic.p_ctg.gfastats

cd ..

## Purge_dups. In order to perform purge_dups there are many steps that have to be done before.

mkdir purgedups
cd purgedups
prim_asm='../hifiasm_HiC/${ASM_NAME}_asm_prim.hic.p_ctg.fasta'
sec_asm='../hifiasm_HiC/${ASM_NAME}_asm_prim.hic.a_ctg.fasta'
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
busco -i purged.fa -o quality_stats/busco --lineage $lineage -c ${THREADS} -m geno

conda deactivate

quast.py hap.fa --large --est-ref-size $genome_size -o quality_stats/quast

echo "Step 2 has been performed. Assembly and its metrics before and after purging the duplicates aswell, check folders assembly_metrics/quality_stats and purgedups/quality_stats"

