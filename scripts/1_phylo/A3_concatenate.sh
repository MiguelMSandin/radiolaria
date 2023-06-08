#!/bin/bash

FILE18S="seqs_18S_align-auto_trimed05.fasta"
FILE28S="seqs_28S_align-auto_trimed05.fasta"
OUT="seqs_concat_align-auto_trimed05.fasta"

fastaConcat.py -f $FILE18S $FILE28S -o $OUT -a
