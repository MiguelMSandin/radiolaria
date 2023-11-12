#!/bin/bash

[ ! -d "data/phylo/" ] && mkdir -p "data/phylo/"
[ ! -d "tmp" ] && mkdir -p "tmp"

# Extract 18S sequences ----------------------------------------------------------------------------
cut -f 2 "seqs_id_raw.tsv" > "tmp/acnu18s"
sed -i 's/acnu_18S//g' "tmp/acnu18s"
sed -i '/^NA$/d' "tmp/acnu18s"
sed -i '/^$/d' "tmp/acnu18s"

# Download
ACNUs=$(tr '\n' ',' < tmp/acnu18s | sed 's/,$//g')
efetch -db nucleotide -id $ACNUs -format fasta > "data/phylo/seqs_18S.fasta"

# Change names
sed -i 's/\..*//g' data/phylo/seqs_18S.fasta

while read ID acnu18s acnu28s taxo; do
	if [ $acnu18s != "NA" ]; then sed -i "s/$acnu18s/$ID-$taxo/g" data/phylo/seqs_18S.fasta; fi
done < "seqs_id_raw.tsv"


# Extract 28S sequences ----------------------------------------------------------------------------
cut -f 3 "seqs_id_raw.tsv" > "tmp/acnu28s"
sed -i 's/acnu_28S//g' "tmp/acnu28s"
sed -i '/^NA$/d' "tmp/acnu28s"
sed -i '/^$/d' "tmp/acnu28s"

# Download
ACNUs=$(tr '\n' ',' < tmp/acnu28s | sed 's/,$//g')
efetch -db nucleotide -id $ACNUs -format fasta > "data/phylo/seqs_28S.fasta"

# Change names
sed -i 's/\..*//g' data/phylo/seqs_28S.fasta

while read ID acnu18s acnu28s taxo; do
	if [ $acnu28s != "NA" ]; then sed -i "s/$acnu28s/$ID-$taxo/g" data/phylo/seqs_28S.fasta; fi
done < "seqs_id_raw.tsv"


rm -rf tmp

