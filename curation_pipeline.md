# Curation of Radiolaria sequences

Miguel M. Sandin  
02-07-2020  
miguelmendezsandin@gmail.com

## Dependencies
- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) 
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) ([Katoh and Standley, 2013](https://academic.oup.com/mbe/article/30/4/772/1073398)).
- [Mothur](https://github.com/mothur/mothur/releases/tag/v.1.44.1) ([Schloss et al., 2009](https://journals.asm.org/doi/10.1128/AEM.01541-09)).
- [python3](https://docs.python.org/3/reference/).
- [R](https://www.r-project.org/).
- [RaxML](https://cme.h-its.org/exelixis/web/software/raxml/index.html) ([Stamatakis, 2014](https://academic.oup.com/bioinformatics/article/30/9/1312/238053)).
- [trimAL](http://trimal.cgenomics.org/) ([Capella-Gutiérrez et al., 2009](https://academic.oup.com/bioinformatics/article/25/15/1972/213148)).
- [vsearch](https://github.com/torognes/vsearch) ([Rognes et al., 2016](https://peerj.com/articles/2584/)).
### Optional
- [AliView](http://www.ormbunkar.se/aliview/#DOWNLOAD) ([Larsson, A. 2014](https://academic.oup.com/bioinformatics/article/30/22/3276/2391211)).
- [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) ([Rambaut 2016](http://tree.bio.ed.ac.uk/software/figtree/)).
- [Rstudio](https://rstudio.com/products/rstudio/download/).
- [SeaView](http://doua.prabi.fr/software/seaview) ([Gouy et al., 2010](https://academic.oup.com/mbe/article/27/2/221/970247)).
  
---  
## Step 1. Building a detailed reference alignment
I have used as reference phylogenetic frameworks the last morpho-molecular classifications of **Acantharea** ([Decelle et al. 2012](https://www.sciencedirect.com/science/article/abs/pii/S1434461011000988)), **Collodaria** ([Biard et al. 2015](https://www.sciencedirect.com/science/article/pii/S1434461015000231)), **Nassellaria** ([Sandin et al. 2019](https://www.sciencedirect.com/science/article/pii/S143446101830110X)) and **Spumellaria** ([Sandin et al. 2021](https://www.sciencedirect.com/science/article/pii/S1434461021000158)). Following these phylogenetic works I have selected 2 to 4-6 representative sequences of each clade (both morphologically described and environmental) to build a reference alignment (in total 211 sequences), in order to gather the most important and detailed genetic diversity. Such sequences were manually selected due to their high quality and their length, comprising the full 18S rDNA gene (when available). As showed in the reference phylogenetic frameworks, the partial 28S (D1+D2) rDNA gene was also used and concatenated in order to increase phylogenetic signal. 
  
### 1.1. Align sequences
Sequences were aligned with MAFFT using the L-INS-i algorithm (`--localpair`) and 1000 iterative refinement cycles for high accuracy:
```
mafft --localpair --maxiterate 1000 FILE > FILE_aligned
```
  
### 1.2. Check alignment
Alignment was manually check and corrected in SeaView (version 4.6.1) in order to align analogous positions and check for possible missalignments. 
  
### 1.3. Trimm alignment
Corrected alignment matrix was trimmed with trimal selecting a 30% threshold:
```
trimal -in FILE_alignedC -out FILE_alignedC_trimmed -gt 0.3
```
  
### 1.4. Phylogenetic assessment
Two (semi-replicate) analysis were used to look for agreement -or not- and resolve possible dubious or low statistically supported groups/patterns using RAxML:
#### 1.4.1 GTR+CAT as substitution model and 1000 bootstraps (BS)
```
raxmlHPC-PTHREADS-SSE3 -T 5 -m GTRCAT -c 25 -p $RANDOM -x $(date +%s) -d -f a -N 1000 -n FILE_alignedC_trimmed_CAT -s FILE_alignedC_trimmed
```
#### 1.4.2. GTR+GAMMA with 1000 BS
```
raxmlHPC-PTHREADS-SSE3 -T 5 -m GTRGAMMA -p $RANDOM -x $(date +%s) -d -f a -N $1000 -n FILE_alignedC_trimmed_GAMMA -s FILE_alignedC_trimmed
```
  
### 1.5. Re-ordering alignment
The resulting alignment was ordered following the phylogenetic tree in order to ease the manual correction of misalignments with the function *[fastaReorder.py](https://github.com/MiguelMSandin/random/blob/main/phylogenetics/fastaReorder.py)* from the repository *[random](https://github.com/MiguelMSandin/random)*.  
Steps **1.2** to **1.5** were repeated until the phylogenetic tree had a consistent topology and in agreement with previous studies ([Decelle et al. 2012](https://www.sciencedirect.com/science/article/abs/pii/S1434461011000988); [Biard et al. 2015](https://www.sciencedirect.com/science/article/pii/S1434461015000231); [Sandin et al. 2019](https://www.sciencedirect.com/science/article/pii/S143446101830110X); [2021](https://www.sciencedirect.com/science/article/pii/S1434461021000158)). Phylogenetic trees were visualized in FigTree (version 1.4.3). Resulting alignment matrix was set as reference alignment for further steps.

<sub>*Note*: Steps **1.1** to **1.3** were carried out independently for 18S and 28S rDNA genes. Then both genes were concatenated for step **1.4** and further.</sub>
  
---
## Step 2. Creating a reference backbone phylogeny
I continued adding the rest of publicly available and morphologically identified sequences detailed in the reference phylogenetic frameworks to the previous reference alignment (from **step 1**) using the function `--add` from MAFFT and the FFT-NS-2 algorithm with a high gap opening penalty (`--op 5`) as follows:
```
mafft --thread 2 --inputorder --op 5.0 --add FILE2 --6merpair --maxiterate 1000 FILE_alignedC > FILE2_aligned
```
And repeating steps **1.2** to **1.5** until consistent phylogenetic tree.
  
---
## Step 3. Download all publicly available environmental sequences
Reference 18S rDNA gene sequences from **step 2** were ‘blasted’ against NCBI to retrieve all publicly available environmental sequences of the 18S rDNA gene that belong to Radiolaria, as detailed in ‘Step 6’ of [EukRef](http://eukref.org/curation-pipeline-overview/). Previous versions of [PR2](https://github.com/pr2database) ([Guillou et al., 2013](https://academic.oup.com/nar/article/41/D1/D597/1064851)) were also blasted against NCBI in order to retrieve sequences not belonging to either Acantharea, Collodaria, Nassellaria or Spumellaria (i.e.; Rad-A, Rad-B, Rad-C or Radiolaria_X). Final dataset was checked in order to remove duplicates and other regions of the rDNA (except for the partial 28S rDNA gene that was left and concatenated, when available, to increase phylogenetic signal).

---
## Step 4. Check for chimeras
Due to the difficulties in chimera detection, I have used two different algorithms implemented in:
  
### 4.1. mothur
```
mothur "#chimera.uchime(fasta=NEW, reference=REFERENCE)"
```
  
### 4.2. and vsearch
```
vsearch --uchime_ref NEW --db REFERENCE --nonchimeras NEW_noChim --borderline NEW_poChim --chimeras NEW_Chim
```
  
Using as REFERENCE sequences the complete PR2 v4.11.0 database and as NEW those sequences extracted from **step 3**. Sequences detected as chimeras by any of the two methods were automatically removed and no longer considered in downstream analysis.
  
  
### 4.3. Blasting
A third in-house method was also used in order to detect as many chimeric or problematic sequences as possible. Here I blasted against NCBI independently the first and the last 300 bp of each sequence and compared the results. If there were less than 20 exact matches among the first 100 matches, the sequence was considered as *maybe* chimeric. If there are no exact matches among the first 100 matches, or it matches only with itself, the sequence was detected as *chimeric*. For further details see the script ‘[chimera_blast.R](https://github.com/MiguelMSandin/Radiolaria_reference/blob/master/scripts/chimera_blast.R)’ (*recommended to open it in Rstudio*). Sequences detected either as maybe chimeric or chimeric were noted and manually checked in downstream phylogenetic assessments.
  
---
## Step 5. Annotating sequences
All environmental sequences resulting from **step 4** were automatically annotated with vsearch taking as reference sequences those morphologically described:
```
vsearch --usearch_global ENVIRONMENTAL --db REFERENCE --blast6out OUTPUT --log OUTPUT_log
```
  
---
## Step 6. Adding environmental sequences to the backbone phylogeny
Environmental annotated sequences from **step 5** were added to the reference alignment created in **step 2**. These sequences were gradually added according to their ‘reliability’ so it is easier to identify chimeric, *dubious* or bad quality sequences: 
1. Firstly, were added those environmental sequences retrieved in previous studies and phylogenetically placed in a reference phylogenetic framework ([Decelle et al. 2012](https://www.sciencedirect.com/science/article/abs/pii/S1434461011000988); [Biard et al. 2015](https://www.sciencedirect.com/science/article/pii/S1434461015000231); [Sandin et al. 2019](https://www.sciencedirect.com/science/article/pii/S143446101830110X); [2021](https://www.sciencedirect.com/science/article/pii/S1434461021000158)). 
2. Secondly, were added these sequences that do not correspond within Acantharea, Collodaria, Nassellaria or Spumellaria (i.e.; Rad-A, Rad-B, Rad-C and Radiolaria_X). 
3. Thirdly, and lastly, the rest of the sequences retrieved.  
  
<sub>\*I considered a dubious sequence a sequence that is clearly different (either in pair-wise similarity or phylogenetically) from the others and with a lack of similar sequences. In other words, if a sequence is different from others but has 2 or (ideally) more similar sequences, and preferably from different studies, is not considered as dubious.</sub>
  
For every group of sequences that were added, the following steps were carried out:
  
### 6.1. Align sequences against previously created reference alignment
```
mafft --thread 2 --inputorder --op 5.0 --add FILEX --6merpair --maxiterate 1000 FILE_alignedC > FILEX_aligned
```
  
### 6.2. Trimming alignment
```
trimal -in FILEX_aligned -out FILEX_aligned_trimmed -gt 0.3
```
  
### 6.3. Phylogenetic assessment
Explore phylogenetic patterns, long branches, branches in unresolved positions, contrasting topology of previously resolved clades, … 
  
#### 6.3.1. RAxML GTR+CAT 100 BS
```
raxmlHPC-PTHREADS-SSE3 -T 6 -m GTRCAT -c 25 -p $RANDOM -x $(date +%s) -d -f a -N 100 -n FILEX_aligned_trimmed_CAT -s  FILEX_aligned_trimmed
```
  
#### 6.3.2. RAxML GTR+GAMMA 100 BS
```
raxmlHPC-PTHREADS-SSE3 -T 6 -m GTRGAMMA -p $RANDOM -x $(date +%s) -d -f a -N $100 -n FILEX_aligned_trimmed_GAMMA -s  FILEX_aligned_trimmed
```
  
### 6.4. Reorder fasta 
Reordering fasta file with the function ‘reorderFastaTree’ from ‘fastaFunctions.R’, taking as reference tree that of 6.3.1 (or 6.3.2).
  
### 6.5. Manual removing of problematic sequences
Manual and final identification of chimeric, dubious or bad quality sequences for its deletion in AliView (version 1.26).
Repeat steps **6.3** to **6.5** until consistent phylogeny. 
  
<sub>*Note*: in order to ease and speed **step 6**, after removing problematic sequences (**6.5**), I directly run the phylogenetic analysis (**6.3**), without realigning the raw sequences (**6.1**) and trimming the alignment (**6.2**). This is why I added a last confirmatory check:</sub>
  
---
## Step 7. Last confirmatory check
All environmental sequences added in **step 6** that passed the stringent criteria were (re-)aligned to the backbone phylogenetic framework from **step 2** for a last phylogenetic assessment as described in steps **6.1** to **6.3**. And, if needed, step **6.3** was repeated after step **6.4** and manual checking and correction of the alignment in AliView (version 1.26).  
  
---
## Step 8. Final taxonomic annotation and corrections
Final phylogenetic trees were used for the manual curation and correction of the taxonomic annotation of the enviromental sequences according to bootstrap support and phylogenetic relatedness to morphologically described sequences. Annotation was done according to the original description of:  
- **Acantharea** ([Decelle et al. 2012](https://www.sciencedirect.com/science/article/abs/pii/S1434461011000988)): Clade A, Clade B, Clade C, Clade D, Clade E, Clade F and the environmental clades Acantharea 1, Acantharea 2 (subclades “a” and “b”), Acantharea 3 (subclades “a” and “b”) and Acantharea 4 (subclades “a” to “d”).  
<sub>*Note*: Environmental clades of Acantharea were changed from roman numbers (I, II, III, IV) to Arabic numbers (1, 2, 3, 4 respectively) in order to homogenize the annotation.</sub>  
- **Collodaria** ([Biard et al. 2015](https://www.sciencedirect.com/science/article/pii/S1434461015000231)): Collophidiidae, Collosphaeridae and Sphaerozoidae, and some environmental sequences (annotated as Collodaria_X) as originally curated by Tristan Biard.  
- **Nassellaria** ([Sandin et al. 2019](https://www.sciencedirect.com/science/article/pii/S143446101830110X)): Acanthodesmoidea (gathering Acanthodesmidae and Triospyrididae), Archipiliidae, Artostrobiidae, Carpocaniidae, Cycladophoridae, Eucyrtidiidae, Lychnocaniidae, Orosphaeridae (as described in [Nakamura et al., 2020](https://www.sciencedirect.com/science/article/pii/S1434461019300707)), Plagiacanthoidea (gathering Cannobotryidae, Lophophaenidae, Sethoperidae and Sethophormidae), Pterocorythidae and Theopiliidae and the environmental clades Nassellaria_1, Nassellaria_2, Nassellaria_3 and Nassellaria_4.  
<sub>*Note1*: Acropyramioidea has no morphologically described representatives for the 18S rDNA gene, yet Nassellaria_1 is phylogenetically related to this clade.</sub>  
<sub>*Note2*: Nassellaria_3 shows very strong phylogenetic relationships with Artostrobiidae.</sub>  
- **Spumellaria** ([Sandin et al. 2021](https://www.sciencedirect.com/science/article/pii/S1434461021000158)): Actinommidae, Centrocuboidea (gathering Centrocubidae and Excentroconchidae), Hexastyloidea (gathering Hexalonchidae, Hexastylidae and Hollandosphaeridae), Spongodiscoidea (gathering Euchitoniidae and Spongobrachiidae), Spongopyloidea (gathering Cristallosphaeridae and Spongopylidae), Pylonioidea (gathering Pyloniidae and Tholoniidae), Astrosphaeridae, Lithelidae, Artiscinae, Rhizosphaeridae, Spongosphaeridae and Stylodictyidae and the envrionmental clades Spumellaria_1, Spumellaria_2, Spumellaria_3, Spumellaria_4 (subclades “a” and “b”) and Spumellaria_6.  
<sub>*Note*: Spumellaria_5 (from [Sandin et al. 2021](https://www.sciencedirect.com/science/article/pii/S1434461021000158)) was emended within Liosphaeroidea (gathering Astrosphaeridae) as “Liosphaeroidea_X”.</sub>  
  
Environmental clades (Rad-A, Rad-B, Rad-C and Radiolaria_X) were left as previously annotated, yet a more detailed taxonomy was applied when clear phylogenetic patterns (bootstrap ≥ 90 and consensus between GTR+GAMMA and GTR+CAT):
- **RAD-A**: Clades from “a” to “g”.  
- **RAD-B**: Group-I (subclades “a” to “f”), group-II, group-III and group-IV (with Taxopodida and subclades “a” to “e” and some nonassignable group of sequences “X”).  
- **RAD-C**: Clades “a” and “b”.  
  
---
## Summary and concluding remarks
In total, 4569 publicly available sequences associated to Radiolaria have been taxonomically curated and annotated [(files/pr2_Radiolaria_MMS.xlsx)](https://github.com/MiguelMSandin/Radiolaria_reference/blob/master/files/pr2_Radiolaria_MMS.xlsx). From these, 4556 are longer than 500bp and have less than 20 ambiguities [(files/pr2_rads_updated.fasta)](https://github.com/MiguelMSandin/Radiolaria_reference/blob/master/files/pr2_rads_updated.fasta) and 289 sequences are new regarding older PR2 versions (as of v4.11.0). Main changes are:  
- **Removing of chimeric, bad quality or dubious sequences**: In total 391 sequences were manually removed, mainly from environmental clades (such as Acantharea_1, Acantharea_3, Rad-B, Radiolaria_X) and Liosphaeroidea (Spumellaria), among others.  
- **Updating Acantharea** sequences after Decelle et al. (2012): 457 sequences.  
- Formal **inclusion of Collodaria within Nassellaria** (after [Nakamura et al., 2020](https://www.sciencedirect.com/science/article/pii/S1434461019300707) and unpublished results from Miguel M. Sandin PhD thesis, 2019).  
- **Updating Nassellaria** sequences after Sandin et al. (2019): 459 sequences, of which 190 formerly belonged to Collodaria.  
- **Updating Spumellaria** sequences after Sandin et al. (2020): 3017 sequences.  
- **Finer taxonomic resolution of environmental clades** RAD-A (107 sequences), RAD-B (455 sequences) and RAD-C (50 sequences).  
- Up to 24 sequences, phylogenetically related to Acantharia, were not possible to be assigned to any previously described clade, and were annotated as “**Radiolaria_X**”.  

Here we built a backbone phylogenetic framework of all the described diversity of Radiolaria to which we have aligned the rest of publicly available environmental sequences, after removing problematic sequences (chimeras, dubious or bad quality). All this gathers fully reliable sequences representing the biggest part of the diversity of Radiolaria, that we can trust and use for phylogenetic reconstruction, phylogenetic placement of short reads or annotation of the still to come full length environmental rDNA sequences from Oxford Nanopore Technologies or Pac-Bio.  
  
These sequences have been deposited in PR2 v5.0.0 (to be published).  
  
Lastly, I would like to mention that **taxonomic hierarchy may differ from one reference to another**. For example in Cavalier-smith ([2017](https://link.springer.com/article/10.1007/s00709-018-1241-1)) Polycystinea (gathering Spumellaria, Nassellaria and Collodaria)  holds a “superclass” position along with Spasmaria (gathering Acantharea and Sticholonche), and  in Adl et al. ([2019](https://onlinelibrary.wiley.com/doi/epdf/10.1111/jeu.12691)) Polycystinea holds a “class” position along with Acantharea and Taxopodida. You may therefore select specific taxa according to the scope of your study bearing in mind these differences in taxonomic hierarchy. 
  
  
---
## References
- Adl SM, Bass D, Lane CE, Lukeš J, Schoch CL, Smirnov A, Agatha S, Berney C, Brown MW, Burki F, Cárdenas P, Čepička I, Chistyakova L, Del Campo J, Dunthorn M, Edvardsen B, Eglit Y, Guillou L, Hampl V, Heiss AA, Hoppenrath M, James TY, Karnkowska A, Karpov S, Kim E, Kolisko M, Kudryavtsev A, Lahr DJG, Lara E, Le Gall L, Lynn DH, Mann DG, Massana R, Mitchell EAD, Morrow C, Park JS, Pawlowski JW, Powell MJ, Richter DJ, Rueckert S, Shadwick L, Shimano S, Spiegel FW, Torruella G, Youssef N, Zlatogursky V, Zhang Q (**2019**) Revisions to the Classification, Nomenclature, and Diversity of Eukaryotes. *J Eukaryot Microbiol*. 66(1):4-119. doi: [10.1111/jeu.12691](https://onlinelibrary.wiley.com/doi/epdf/10.1111/jeu.12691)  
- Biard T, Bigeard E, Audic S, Poulain J, Gutierrez-Rodriguez A, Pesant S, Stemmann L, Not F (**2017**) Biogeography and diversity of Collodaria (Radiolaria) in the global ocean. *ISME J*. 11(6):1331-1344. doi: [10.1038/ismej.2017.12](https://www.nature.com/articles/ismej201712)  
- Biard T, Pillet L, Decelle J, Poirier C, Suzuki N, Not F. Towards an Integrative Morpho-molecular Classification of the Collodaria (Polycystinea, Radiolaria). Protist. 2015 Jul;166(3):374-88. doi: [10.1016/j.protis.2015.05.002](https://www.sciencedirect.com/science/article/pii/S1434461015000231)
- Cavalier-Smith T, Chao EE, Lewis R (**2018**) Multigene phylogeny and cell evolution of chromist infrakingdom Rhizaria: contrasting cell organisation of sister phyla Cercozoa and Retaria. *Protoplasma*. 255(5):1517-1574. doi: [10.1007/s00709-018-1241-1](https://link.springer.com/article/10.1007/s00709-018-1241-1)  
- Capella-Gutiérrez S, Silla-Martínez JM, Gabaldón T (**2009**) trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. *Bioinformatics*. 25(15):1972-3. doi: [10.1093/bioinformatics/btp348](https://academic.oup.com/bioinformatics/article/25/15/1972/213148)  
- Decelle J, Suzuki N, Mahé F, Vargas Cde, Not F (**2012b**) Molecular Phylogeny and Morphological Evolution of the Acantharea (Radiolaria). *Protist*. 163, 435–450. doi: [10.1016/j.protis.2011.10.002](https://www.sciencedirect.com/science/article/abs/pii/S1434461011000988)  
- Gouy M, Guindon S, Gascuel O (**2010**) SeaView Version 4: A Multiplatform Graphical User Interface for Sequence Alignment and Phylogenetic Tree Building. *Mol. Biol. Evol*. 27, 221–224. doi: [10.1093/molbev/msp259](https://academic.oup.com/mbe/article/27/2/221/970247)  
- Guillou L, Bachar D, Audic S, Bass D, Berney C, Bittner L, Boutte C, Burgaud G, de Vargas C, Decelle J, Del Campo J, Dolan JR, Dunthorn M, Edvardsen B, Holzmann M, Kooistra WH, Lara E, Le Bescot N, Logares R, Mahé F, Massana R, Montresor M, Morard R, Not F, Pawlowski J, Probert I, Sauvadet AL, Siano R, Stoeck T, Vaulot D, Zimmermann P, Christen R (**2013**) The Protist Ribosomal Reference database (PR2): a catalog of unicellular eukaryote small sub-unit rRNA sequences with curated taxonomy. *Nucleic Acids Res*. 41(Database issue):D597-604. doi: [10.1093/nar/gks1160](https://academic.oup.com/nar/article/41/D1/D597/1064851)  
- Katoh K, Standley DM (**2013**) MAFFT multiple sequence alignment software version 7: Improvements in performance and usability. *Mol. Biol. Evol*. 30, 772–780. doi: [10.1093/molbev/mst010](https://academic.oup.com/mbe/article/30/4/772/1073398)  
- Larsson A (**2014**) AliView: a fast and lightweight alignment viewer and editor for large data sets. *Bioinformatics*. 30(22): 3276-3278. doi: [10.1093/bioinformatics/btu531](https://academic.oup.com/bioinformatics/article/30/22/3276/2391211)  
- Nakamura Y, Sandin MM, Suzuki N, Somiya R, Tuji A, Not F (**2020**) Phylogenetic revision of the order Entactinaria - Paleozoic relict Radiolaria (Rhizaria, SAR). *Protist* 171(1):125712. doi: [10.1016/j.protis.2019.125712](https://www.sciencedirect.com/science/article/pii/S1434461019300707)  
- R Core Team (**2019**). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL [https://www.R-project.org/](https://www.R-project.org/)  
- Rambaut A (**2016**) FigTree version 1.4.3. [http://tree.bio.ed.ac.uk/software/figtree/](http://tree.bio.ed.ac.uk/software/figtree/)  
- Rognes T, Flouri T, Nichols B, Quince C, Mahé F (**2016**) VSEARCH: a versatile open source tool for metagenomics. *PeerJ* 4, e2584. doi: [10.7717/peerj.2584](https://peerj.com/articles/2584/)  
- Sandin MM, Biard T, Romac S, O'Dogherty L, Suzuki N, Not F (**2021**) A Morpho-molecular Perspective on the Diversity and Evolution of Spumellaria (Radiolaria). *Protist*. 172(3):125806. doi: [10.1016/j.protis.2021.125806](https://www.sciencedirect.com/science/article/pii/S1434461021000158)   
- Sandin MM, Pillet L, Biard T, Poirier C, Bigeard E, Romac S, Suzuki N, Not F (**2019**) Time Calibrated Morpho-molecular Classification of Nassellaria (Radiolaria). *Protist* 170, 187–208. doi: [10.1016/j.protis.2019.02.002](https://www.sciencedirect.com/science/article/pii/S143446101830110X)  
- Schloss PD, Westcott SL, Ryabin T, Hall JR, Hartmann M, Hollister EB, Lesniewski RA, Oakley BB, Parks DH, Robinson CJ, Sahl JW, Stres B, Thallinger GG, Van Horn DJ, Weber CF (**2009**) Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. *Appl. Environ. Microbiol*. 75, 7537–7541. doi:[10.1128/AEM.01541-09](https://journals.asm.org/doi/10.1128/AEM.01541-09)  
- Stamatakis A (**2014**) RAxML version 8: A tool for phylogenetic analysis and post-analysis of large phylogenies. *Bioinformatics* 30, 1312–1313. doi: [10.1093/bioinformatics/btu033](https://academic.oup.com/bioinformatics/article/30/9/1312/238053)   
