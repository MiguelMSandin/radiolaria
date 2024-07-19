#!/bin/bash
#
#SBATCH -o slurm.%x.%N.%j.out
#SBATCH -e slurm.%x.%N.%j.err
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user YOUR_EMAIL_ADDRESS
#SBATCH --partition fast
#SBATCH --mem 30GB
#SBATCH -t 1-0:00
#SBATCH -J mtMC2vl

module load paml/4.9

CONTROL="all_filtered_align-linsi_trim05_MC10_step2.ctl"

mcmctree_NS2000 $CONTROL

# mv "SeedUsed" ${CONTROL/.ctl/.seed}
# mv "FigTree.tre" ${CONTROL/.ctl/_mcmc.tre}
mv "FigTree.tre" "all_filtered_align-linsi_trim05_2311_mcmcTree_MC10_step1vl_b0.tre"
