          seed = -1
       seqfile = all_filtered_align-linsi_trim05_2311.fasta
      treefile = all_filtered_align-linsi_trim05_2311_raxmlCAT_mcmcTree_MC10_calibrated.tre
      mcmcfile = all_filtered_align-linsi_trim05_2311_mcmcTree_MC10_step1.txt
       outfile = all_filtered_align-linsi_trim05_2311_mcmcTree_MC10_step1.log

         ndata = 1
       seqtype = 0    * 0: nucleotides; 1:codons; 2:AAs
       usedata = 2    * 0: no data; 1:seq like; 2:normal approximation; 3:out.BV (in.BV)
         clock = 2    * 1: global clock; 2: independent rates; 3: correlated rates
       RootAge = '>11<22'  * safe constraint on root age, used if no fossil for root.

       BDparas = 1 1 0  * birth, death, sampling
   kappa_gamma = 6 2      * gamma prior for kappa
   alpha_gamma = 1 1      * gamma prior for alpha

   rgene_gamma = 2 40 0 1 * gammaDir prior for rate for genes
  sigma2_gamma = 1 10 2   * gammaDir prior for sigma^2     (for clock=2 or 3)

         print = 1   * 0: no mcmc sample; 1: everything except branch rates 2: everything; -1: for already run mcmc sample
        burnin = 250000
      sampfreq = 2
       nsample = 500000

    checkpoint = 1 * 0: nothing; 1 : save; 2: resume

*** Note: Make your window wider (100 columns) before running the program.
