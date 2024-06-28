# Pipeline used for Whole Genome Assembly
Proposed pipeline for de novo assembly. This first version is prepared for using HIFI and HIC reads.

# 1. Genome size,  heterozygosity and ploidy estimation (1_reads_assess.sh)
Tools used:
* NanoPlot, installed through git clone: https://github.com/wdecoster/NanoPlot
* LongQC, installed through git clone: https://github.com/yfukasawa/LongQC
* kmc, installed through apt: https://github.com/refresh-bio/KMC
* GenomeScope2. Installed using mamba: https://github.com/tbenavi1/genomescope2.0
* FastK installed through git clone: https://github.com/thegenemyers/FASTK
* Smudgeplot. Installed through git clone: https://github.com/KamilSJaron/smudgeplot/tree/sploidyplot
      
# 2. Assembly quality check and purge duplicates (2_asm_assess.sh) 
Tools used:
* hifiasm, installed through git clone:  https://github.com/chhylp123/hifiasm
* BUSCO, installed via mamba.  mamba install bioconda::busco
* quast.py, installed  through git clone (check requirements!)  https://github.com/ablab/quast.git
* gfastats, installed through git clone: https://github.com/vgl-hub/gfastats
* purge_dups, installed through git clone: https://github.com/dfguan/purge_dups
    
# 3. First decontamination (3_pre_scf_decon.sh)
Tools used:
* FCS, https://github.com/ncbi/fcs
* Whokaryote, https://github.com/LottePronk/whokaryote
      
# 4. Scaffolding (4_scf.sh)
Tools used:
* Yahs, installed through git clone: https://github.com/c-zhou/yahs
* bwa, installed through apt
* samtools, installed through apt
* two_read_bam_combiner.pl, from Arima Genomics Pipeline,  https://github.com/ArimaGenomics/mapping_pipeline

# 5. Second decontamination (5_blob_post_scf.sh)
Tools used:
* Blobtoolkit
  
# 6. Manual curation (6_pret_post_scf.sh)
Tools used:
* PretextMap
* PretextView
* PretextGraph
  
# 7. Post Pretext processing (7_post_pret.sh)

    
