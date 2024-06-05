##2 ASSEMBLY STEP AND ASSEMBLY QUALITY CHECK
## usage ./assembly_quality_check.sh "HIFI_reads" "HIC_1_reads" "HIC_2_reads" "cpu"
HIFI = ${1}
HIC1 = ${2}
HIC2 = ${3}
THR = ${4}
mkdir hifiasm_HiC
cd hifiasm_HiC
hifiasm -o assembly_hic_prim -t ${THR} --h1 ${HIC1} --h2 ${HIC2} ${HIFI} --primary

# get fasta files from hifiasm assemblies
awk '/^S/{print ">"$2;print $3}' assembly_hic_prim.hic.p_ctg.gfa > asm_prim.hic.p_ctg.fasta
awk '/^S/{print ">"$2;print $3}' assembly_hic_prim.hic.a_ctg.gfa > asm_prim.hic.a_ctg.fasta
cd ..

## Metrics using QUAST + Busco + gfastats

# 3.1 BUSCO of contigs obtained from assembler
mkdir -p assembly_metrics/busco_out
cd assembly_metrics

#mamba activate busco
busco -i ../asm_prim.hic.p_ctg.fasta -o busco_out --lineage arthropoda_odb10 -c ${THR} -m geno

#mamba deactivate
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

mamba activate busco
# Applying BUSCO to the first 3 steps output.
busco -i purged.fa -o quality_stats/busco --lineage arthropoda_odb10 -c ${THR} -m geno

mamba deactivate

quast.py hap.fa --large --est-ref-size 1180556766 -o quality_stats/quast


