#!/bin/bash

FILE="seqs_concat_align-auto_trimed05.fasta"
PREFIX="seqs_concat_align-auto_trimed05.iqtree_gtrGi"

BS="100"
MEM="10G"
RUNS="10"
THREADS="2"

iqtree -s $FILE -st DNA -pre $PREFIX -b $BS -seed $(date +%s) -mem $MEM -nt $THREADS --runs $RUNS -wbtl
