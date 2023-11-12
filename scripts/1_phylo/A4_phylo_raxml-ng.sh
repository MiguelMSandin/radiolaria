#!/bin/bash

FILE=""
PREFIX=""

raxml-ng --parse --msa $FILE --model GTR+G --prefix $PREFIX.parser

BS="1000"
THREADS=""
START_TREES="100"

raxml-ng --all --msa $PREFIX.parser.raxml.rba --tree pars{$START_TREES} --prefix $PREFIX --seed $RANDOM --threads $THREADS --bs-trees $BS
