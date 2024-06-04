##2 ASSEMBLY STEP AND ASSEMBLY QUALITY CHECK

HIFI = ${1}
HIC1 = ${2}
HIC2 = ${3}
THR = ${4}
mkdir hifiasm_HiC
cd hifiasm_HiC
hifiasm -o assembly_hic_prim -t ${4}--h1 ${2} --h2 ${3} ${1} --primary

# get fasta files from hifiasm assemblies
awk '/^S/{print ">"$2;print $3}' assembly_hic_prim.hic.p_ctg.gfa > asm_prim.hic.p_ctg.fasta
awk '/^S/{print ">"$2;print $3}' assembly_hic_prim.hic.a_ctg.gfa > asm_prim.hic.a_ctg.fasta
cd ..

## Step 3. Metrics using QUAST + Busco

# 3.1 BUSCO of contigs obtained from the different assemblers
mkdir -p assembly_metrics/busco_out
cd assembly_metrics

#mamba activate busco
busco -i ../asm_prim.hic.p_ctg.fasta -o busco_out --lineage arthropoda_odb10 -c ${4} -m geno

#mamba deactivate
# 3.2 QUAST metrics of the purge dups 3rd step.
mkdir -p quast
quast.py ../asm_prim.hic.p_ctg.fasta --large --est-ref-size 1180556766 -o quast

/opt/gfastats/build/bin/gfastats asm_prim.hic.p_ctg.fasta > asm_prim.hic.p_ctg.gfastats


