#!/bin/bash
#
#SBATCH -o slurm.%x_%a.%N.%j.out
#SBATCH -e slurm.%x_%a.%N.%j.err
#SBATCH --mail-type ALL
#SBATCH --mail-user YOUR_EMAIL_ADDRESS
#SBATCH --partition long
#SBATCH --mem 6GB
#SBATCH -t 20-0:00
#SBATCH -J b2MC10i
#SBATCH --array=1-4

module load beast2/2.6.3

FILES=$(ls | grep "xml$")
FILE=$(echo $FILES | cut -d ' ' -f $SLURM_ARRAY_TASK_ID)

beast -overwrite -seed $RANDOM -beagle_SSE $FILE

# THREADS="5"
# beast -threads $THREADS -seed $(date +%s) -beagle_SSE $FILE

#######################################################################
# Make sure the number of arrays fits with the number of *.xml" files!!
# ls | grep -c "xml$"
#######################################################################

