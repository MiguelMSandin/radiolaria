# Beyond the stars of the ocean: Diversity and evolution of extant Radiolaria
  
In this repository you will find all methods, resources and scripts used and described in the following paper:
**REFERENCE**  
  
Briefly:  
## Taxonomic curation of environmental sequences associated to Radiolaria  
All environmental 18S rDNA sequences publicly available (as July 2020) were taxonomically curated as detailed in the [curation_pipeline.md](https://github.com/MiguelMSandin/radiolaria/blob/master/curation_pipeline.md) and publicly accessible in the Protist Ribosomal Reference (PR2) database (from [v4.14.0](https://github.com/pr2database/pr2database/releases/tag/v4.14.0)).  
  
## Alignments for phylogenetic analyses  
In the folder [alignments](https://github.com/MiguelMSandin/radiolaria/tree/master/alignments) there are:  
- A [concatenated alignment](https://github.com/MiguelMSandin/radiolaria/blob/master/alignments/all_filtered_align-linsi_trim05.fasta.gz) used to infer the main phylogeny of Radiolaria (trimmed).  
- A raw alignment containing [18S rDNA](https://github.com/MiguelMSandin/radiolaria/blob/master/alignments/all_18S_filtered_align-linsi.fasta.gz) sequences for future studies. This alignment was obtained with [MAFFT](https://mafft.cbrc.jp/alignment/software/) (using a L-INS-i algorithm and 1000 refinement cycles) and after triming (with [trimAl](http://trimal.cgenomics.org/) and a 5% gap threshold) was used to generate the concatenated alignment.  
- A raw alignment containing [28S rDNA](https://github.com/MiguelMSandin/radiolaria/blob/master/alignments/all_28S_filtered_align-linsi.fasta.gz) sequences for future studies. This alignment was obtained as for the 18S alignment.  
  
## Editable trees resulted from phylogenetic and molecular-clock analyses  
Here there are the phylogenetic trees in nexus format of the different phylogenetic analyses implemented in [RAxML](https://github.com/stamatak/standard-RAxML) under the nucleotide substitution model GTR+CAT over 1000 bootstraps, [RAxML-ng](https://github.com/amkozlov/raxml-ng) under the substitution model GTR+G and a third approach implemented in [IQ-Tree](http://www.iqtree.org/) under the model GTR+F+R10 (chosen based on the highest Bayesian Information Criterion from modelFinder):  
- [all_filtered_align-linsi_trim05_raxmlCAT.tre](https://github.com/MiguelMSandin/radiolaria/tree/master/trees/all_filtered_align-linsi_trim05_raxmlCAT.tre)  
- [all_filtered_align-linsi_trim05_raxml-ngGTR.tre](https://github.com/MiguelMSandin/radiolaria/tree/master/trees/all_filtered_align-linsi_trim05_raxml-ngGTR.tre)  
- [all_filtered_align-linsi_trim05_iqtree.tre](https://github.com/MiguelMSandin/radiolaria/tree/master/trees/all_filtered_align-linsi_trim05_iqtree.tre)
  
There are also available the fossil-calibrated trees obtained with [MCMCTree](http://abacus.gene.ucl.ac.uk/software/MCMCtree.Tutorials.pdf) from the [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) package and with [BEAST2](http://www.beast2.org/):
- [all_filtered_align-linsi_trim05_mc09_mcmcTree.tre](https://github.com/MiguelMSandin/radiolaria/tree/master/trees/all_filtered_align-linsi_trim05_mc09_mcmcTree.tre)
- [all_filtered_align-linsi_trim05_mc09_beast2.tre](https://github.com/MiguelMSandin/radiolaria/tree/master/trees/all_filtered_align-linsi_trim05_mc09_beast2.tre)
  
## Resources  
In this folder there are the [control](https://github.com/MiguelMSandin/radiolaria/tree/master/resources/molClock/mcmcTree) and [xml](https://github.com/MiguelMSandin/radiolaria/blob/master/resources/molClock/beast2/all_filtered_align-linsi_trim05_raxmlCAT_mc09_beast2.xml.gz) files to replicate molecular clock analyses, as well as a tsv version of the [calibration model](https://github.com/MiguelMSandin/radiolaria/blob/master/resources/molClock/fossil_calibrations.tsv), and the [metadata](https://github.com/MiguelMSandin/radiolaria/blob/master/resources/metabarcoding/metadata_assembled_nonRedundant.tsv) for metabarcoding analyses.  
  
## scripts  
And finally this folder contains all scripts used in this study organized by main analyses, as [phylogenetic](https://github.com/MiguelMSandin/radiolaria/tree/master/scripts/1_phylo), [molecular clock](https://github.com/MiguelMSandin/radiolaria/tree/master/scripts/2_molClock) and [metabarcoding](https://github.com/MiguelMSandin/radiolaria/tree/master/scripts/3_metaB) analyses.  
