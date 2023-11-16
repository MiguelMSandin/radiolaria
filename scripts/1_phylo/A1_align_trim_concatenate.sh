#!/bin/bash

FILE18S=""
FILE28S=""
THREADS=""

# Align files --------------------------------------------------------------------------------------
# Using the FFT-NS-i algorithm
# Iterative refinement method
mafft --thread $THREADS --retree 2 --maxiterate 1000 $FILE18S > ${FILE18S/.fasta/_align-fftnsi.fasta}
mafft --thread $THREADS --retree 2 --maxiterate 1000 $FILE28S > ${FILE28S/.fasta/_align-fftnsi.fasta}

# Using the L-INS-i algorithm
# Probably most accurate; recommended for <200 sequences; iterative refinement method incorporating local pairwise alignment information
mafft --thread $THREADS --localpair --maxiterate 1000 $FILE18S > ${FILE18S/.fasta/_align-linsi.fasta}
mafft --thread $THREADS --localpair --maxiterate 1000 $FILE28S > ${FILE28S/.fasta/_align-linsi.fasta}

# trim files ---------------------------------------------------------------------------------------
GAPTHRESHOLD="05"

for FILE in $(ls *align*); do
	trimal -in ${FILE/.fasta/_align-linsi.fasta} -out ${FILE/.fasta/_align-linsi_trim$GAPTHRESHOLD.fasta} -gt 0.$GAPTHRESHOLD
done

# Concatenate --------------------------------------------------------------------------------------
fastaConcat.py -f ${FILE18S/.fasta/_align-fftnsi_trim$GAPTHRESHOLD.fasta} ${FILE28S/.fasta/_align-fftnsi_trim$GAPTHRESHOLD.fasta} -o "all_filtered_align-linsi_trim$GAPTHRESHOLD.fasta" -a
fastaConcat.py -f ${FILE18S/.fasta/_align-linsi_trim$GAPTHRESHOLD.fasta} ${FILE28S/.fasta/_align-linsi_trim$GAPTHRESHOLD.fasta} -o "all_filtered_align-linsi_trim$GAPTHRESHOLD.fasta" -a
