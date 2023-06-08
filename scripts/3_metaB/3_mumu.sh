#!/bin/bash

# Setting names ------------------------------------------------------------------------------------

# Input
FASTA="otus/radiolaria.otus.v20171106.fasta"
OTU_TABLE="otus/radiolaria.otus.v20171106_table.tsv"

# Output
MATCH_LIST="otus/radiolaria.otus.v20171106_similarities.tsv"
NEW_OTU_TABLE="otus/radiolaria_otus_mumu_abun.tsv"
LOG="otus/radiolaria_otus_mumu.log"

# Running vsearch for similarities -----------------------------------------------------------------

vsearch \
	--threads 2 \
	--allpairs_global $FASTA \
	--self \
	--id 0.80 \
	--userfields query+target+id \
	--userout $MATCH_LIST

# awk -F"\t" '$3>0.90' $MATCH_LIST > ${MATCH_LIST/.tsv/_hig.tsv}

# Running mumu for post-clustering -----------------------------------------------------------------

mumu \
	--otu_table $OTU_TABLE \
	--match_list $MATCH_LIST \
	--new_otu_table $NEW_OTU_TABLE \
	--log $LOG

