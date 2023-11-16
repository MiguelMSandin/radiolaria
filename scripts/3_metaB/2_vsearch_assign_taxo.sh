#!/bin/bash

DB="/home/mmsandin/DB/PR2/pr2_version_4.14.0_SSU_taxo_long_euks.fasta"
FASTA="otus/radiolaria.otus.v20171106.fasta"
OUT="otus/radiolaria_otusTaxo_PR2v4.14rads_211216.tsv"

# Extract Radiolaria sequences from DB
sequenceSelect.py -f $DB -p Radiolaria -o ${DB/.fasta/_rads.fasta} -v

# Perform taxonomic assignation
vsearch --thread 2 --usearch_global $FASTA --db ${DB/.fasta/_rads.fasta} --id 0.2 --blast6out $OUT --log ${OUT/.tsv/.log}
