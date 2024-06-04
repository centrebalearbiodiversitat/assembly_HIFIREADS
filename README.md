# assembly
Pipeline used for Whole Genome Assembly

Tools used in genome_metrics.sh

- kmc, installed through apt: https://github.com/refresh-bio/KMC
- GenomeScope2. Installed using mamba: https://github.com/tbenavi1/genomescope2.0
- Smudgeplot. Installed using mamba and the same environment as genomescope2: https://github.com/KamilSJaron/smudgeplot?tab=readme-ov-file

Tools used in assembly_quality_check.sh
  - hifiasm, installed through git clone:  https://github.com/chhylp123/hifiasm
  - BUSCO, installed via mamba.  mamba install bioconda::busco
  - quast.py, installed  through git clone (check requirements!)  https://github.com/ablab/quast.git
  - gfastats, installed through git clone: https://github.com/vgl-hub/gfastats
