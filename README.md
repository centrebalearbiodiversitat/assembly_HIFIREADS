# Pipeline used for Whole Genome Assembly
Proposed pipeline for de novo assembly. This first version is prepared for using HIFI and HIC reads.

# 1. Genome size,  heterozygosity and ploidy estimation (genome_metrics.sh)
   Tools used:
      - NanoPlot, installed through git clone: https://github.com/wdecoster/NanoPlot
      - LongQC, installed through git clone: https://github.com/yfukasawa/LongQC
      - kmc, installed through apt: https://github.com/refresh-bio/KMC
      - GenomeScope2. Installed using mamba: https://github.com/tbenavi1/genomescope2.0
      - Smudgeplot. Installed using mamba and the same environment as genomescope2: https://github.com/KamilSJaron/smudgeplot?tab=readme-ov-file
      
# 2. Assembly quality check and purge duplicates (assembly_quality_check.sh) 
    Tools used:
      - hifiasm, installed through git clone:  https://github.com/chhylp123/hifiasm
      - BUSCO, installed via mamba.  mamba install bioconda::busco
      - quast.py, installed  through git clone (check requirements!)  https://github.com/ablab/quast.git
      - gfastats, installed through git clone: https://github.com/vgl-hub/gfastats
      - purge_dups, installed through git clone: https://github.com/dfguan/purge_dups
    
# 3. First decontamination (decontamination-1.sh)
    Tools used:
      - FCS
      - Whokaryote
      
# 4. Scaffolding (scaffolding.sh)
   Tools used:
    - Yahs, installed through git clone: https://github.com/c-zhou/yahs
    - bwa, installed through apt
    - samtools, installed through apt
    - two_read_bam_combiner.pl, from Arima Genomics Pipeline,  https://github.com/ArimaGenomics/mapping_pipeline

# 5. Second decontamination (decontamination-2.sh)
  Tools used:
  - LRBinner

# 6. Manual curation (getting ready for PretextMap)

    
