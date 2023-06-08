#!/bin/bash

FILE=""
PREFIX=""

BS="100"
MEM="10G"
RUNS="10"
THREADS="2"

iqtree -s $FILE -st DNA -pre $PREFIX -g $CONSTRAINT -m GTR+I+G -b $BS -seed $(date +%s) -mem $MEM -nt $THREADS --runs $RUNS -wbtl -con -net
