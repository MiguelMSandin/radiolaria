#!/bin/bash

# Extarct Radiolaria OTUs
[ ! -d "raw" ] && mkdir -p "raw"
mv "globaldataset.otus.v20171106.tsv" "raw/."
head -n 1 "raw/globaldataset.otus.v20171106.tsv" > "headers"

grep -E "Acantharea|Collodaria|Nassellaria|RAD-|Spumellaria|polycystinea|retaria" "raw/globaldataset.otus.v20171106.tsv" > "tmp"

cat "headers" "tmp" > "raw/radiolaria.otus.v20171106.tsv"

rm -f "headers"

# Extract fasta file
[ ! -d "otus" ] && mkdir -p "otus"
FASTA="otus/radiolaria.otus.v20171106.fasta"

cut -f 1,5 "raw/radiolaria.otus.v20171106.tsv" > "tmp"

sed -i '1d' "tmp"

while read ID SEQ; do
	echo ">$ID" >> $FASTA
	echo "$SEQ" >> $FASTA
done < "tmp"

rm -f "tmp"

# Extract OTU table
cut -f 1,16- "raw/radiolaria.otus.v20171106.tsv" > "otus/radiolaria.otus.v20171106_table.tsv"

# Extract Characteristics table
cut -f 1,4,5,10,11,12,13,14,15 "raw/radiolaria.otus.v20171106.tsv" > "otus/radiolaria.otus.v20171106_chars.tsv"
