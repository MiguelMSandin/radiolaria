#!/bin/bash

FILE18S="all_18S_filtered_align-einsi.fasta"
FILE28S="all_28S_filtered_align-einsi.fasta"
OUT="all_filtered_align-einsi.fasta"

GAPTHRESHOLD="05"

trimal -in $FILE18S -out ${FILE18S/.fasta/_trim$GAPTHRESHOLD.fasta} -gt 0.$GAPTHRESHOLD
trimal -in $FILE28S -out ${FILE28S/.fasta/_trim$GAPTHRESHOLD.fasta} -gt 0.$GAPTHRESHOLD

fastaConcat.py -f $FILE18S $FILE28S -o $OUT -a
fastaConcat.py -f ${FILE18S/.fasta/_trim$GAPTHRESHOLD.fasta} ${FILE28S/.fasta/_trim$GAPTHRESHOLD.fasta} -o ${OUT/.fasta/_trim$GAPTHRESHOLD.fasta} -a

for FILE in $(ls | grep align.*fasta); do
	fastaStats.py -f $FILE >> "fasta_stats.log"
done
