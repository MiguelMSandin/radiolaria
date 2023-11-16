#!/bin/bash

FILE=""
BS="100"
THREADS="2"

let "A=$(echo -n $BS | wc -c) -1"
BS_COUNT="BS10"$A

OUTPUT=${FILE/.fasta/_gtrg_$BS_COUNT}

raxmlHPC-PTHREADS-SSE3 -T $THREADS -m GTRGAMMA -p $RANDOM -x $(date +%s) -d -f a -N $BS -n $OUTPUT -s $FILE

