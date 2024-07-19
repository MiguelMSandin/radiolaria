#!/bin/bash
#
#SBATCH -o slurm.%x.%N.%j.out
#SBATCH -e slurm.%x.%N.%j.err
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user YOUR_EMAIL_ADDRESS
#SBATCH --partition fast
#SBATCH --mem 100GB
#SBATCH -t 1-00:00:00
#SBATCH -J bMC10mg

# module load beast2/2.6.3

PREFIX="all_filtered_align-linsi_trim05_raxmlCAT_BEAST2_MC10"

BURN="25"
OUT="${PREFIX}"
TREE="${OUT}_mcmc.tre"

gunzip *log.gz
# Combine log files --------------------------------------------------------------------------------
/home/umr7144/dipo/mmendezsandin/softwares/beast2/beast/bin/logcombiner \
	-log "$PREFIX-1.log" \
	-log "$PREFIX-2.log" \
	-log "$PREFIX-3.log" \
	-log "$PREFIX-4.log" \
	-b $BURN \
	-o "${OUT}.log"
gzip "${OUT}.log"

gunzip *trees.gz
# Combine tree files -------------------------------------------------------------------------------
/home/umr7144/dipo/mmendezsandin/softwares/beast2/beast/bin/logcombiner \
	-log "$PREFIX-1.trees" \
	-log "$PREFIX-2.trees" \
	-log "$PREFIX-3.trees" \
	-log "$PREFIX-4.trees" \
	-b $BURN \
	-o "${OUT}.trees"

# Run treeannotator --------------------------------------------------------------------------------
/home/umr7144/dipo/mmendezsandin/softwares/beast2/beast/bin/treeannotator -lowMem -heights "median" "${OUT}.trees" $TREE
gzip "${OUT}.trees"
