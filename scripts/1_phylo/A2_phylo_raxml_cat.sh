#!/bin/bash

FILE="seqs_concat_align-ginsi_trimed05.fasta"
BS="100"
THREADS="2"

let "A=$(echo -n $BS | wc -c) -1"
BS_COUNT="BS10"$A

OUTPUT=${FILE/.fasta/_cat_$BS_COUNT}

raxmlHPC-PTHREADS-SSE3 -T $THREADS -m GTRCAT -c 25 -p $RANDOM -x $(date +%s) -d -f a -N $BS -n $OUTPUT -s $FILE

