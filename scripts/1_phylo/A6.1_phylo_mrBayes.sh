
set autoclose=yes nowarnings=yes
execute all_filtered_align-linsi_trim05.nexus
lset nst=6 rates=gamma
mcmc ngen=10000000 Nruns=4 savebrlens=yes file=all_filtered_align-linsi_trim05_mrBayes
sump
sumt
quit
