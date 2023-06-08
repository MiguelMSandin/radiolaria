#!/bin/bash

# 18S file -----------------------------------------------------------------------------------------
FILE="all_subset2_28S.fasta"
GAPTHRESHOLD="05"
THREADS="2"

mafft --thread $THREADS --localpair --maxiterate 1000 $FILE > ${FILE/.fasta/_align-linsi.fasta}

trimal -in ${FILE/.fasta/_align-linsi.fasta} -out ${FILE/.fasta/_align-linsi_trim$GAPTHRESHOLD.fasta} -gt 0.$GAPTHRESHOLD

