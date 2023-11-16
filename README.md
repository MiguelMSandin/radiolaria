# Beyond the stars of the ocean: Diversity and evolution of extant Radiolaria
  
In this repository you will find all methods, resources and scripts used and described in the following paper:
**REFERENCE**  
  
Briefly:  
## Taxonomic curation of environmental sequences associated to Radiolaria  
All environmental 18S rDNA sequences publicly available (as July 2020) were taxonomically curated as detailed in the [curation_pipeline.md](https://github.com/MiguelMSandin/radiolaria/blob/master/curation_pipeline.md) and publicly accessible in the Protist Ribosomal Reference (PR2) database (from [v4.14.0](https://github.com/pr2database/pr2database/releases/tag/v4.14.0)).  
  
## Alignments for phylogenetic analyses  
in the folder [alignments](https://github.com/MiguelMSandin/radiolaria/tree/master/alignments) there are:  
- An alignment containing [18S rDNA](https://github.com/MiguelMSandin/radiolaria/blob/master/alignments/all_18S_filtered_align-linsi.fasta.gz) sequences (untrimmed) for future studies.  
- An alignment containing [28S rDNA](https://github.com/MiguelMSandin/radiolaria/blob/master/alignments/all_28S_filtered_align-linsi.fasta.gz) sequences (untrimmed) for future studies.  
- The [concatenated alignment](https://github.com/MiguelMSandin/radiolaria/blob/master/alignments/all_filtered_align-linsi_trim05.fasta.gz) used to infer the main phylogeny of Radiolaria (trimmed).  
  
## Editable trees resulted from phylogenetic and molecular-clock analyses  
- all_filtered_align-linsi_trim05_raxmlCAT.tre
- all_filtered_align-linsi_trim05_raxml-ngGTR.tre
- all_filtered_align-linsi_trim05_iqtree.tre
- all_filtered_align-linsi_trim05_mc09_mcmcTree.tre
- all_filtered_align-linsi_trim05_mc09_beast2.tre

## Resources  
In this folder there are the [control](https://github.com/MiguelMSandin/radiolaria/tree/master/resources/molClock/mcmcTree) and [xml](https://github.com/MiguelMSandin/radiolaria/blob/master/resources/molClock/beast2/all_filtered_align-linsi_trim05_raxmlCAT_mc09_beast2.xml.gz) files to replicate molecular clock analyses, as well as a editable version of the calibration model, and the [metadata](https://github.com/MiguelMSandin/radiolaria/blob/master/resources/metabarcoding/metadata_assembled_nonRedundant.tsv) for metabarcoding analyses.

## scripts  
And finally this folder contains all scripts used in this study, organized by analyses.
