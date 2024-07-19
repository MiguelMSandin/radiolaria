          seed = -1
       seqfile = all_filtered_align-linsi_trim05_2311.fasta
      treefile = all_filtered_align-linsi_trim05_2311_raxmlCAT_mcmcTree_MC10_calibrated.tre
      mcmcfile = all_filtered_align-linsi_trim05_2311_mcmcTree_MC10_step0.txt
       outfile = all_filtered_align-linsi_trim05_2311_mcmcTree_MC10_step0.log

         ndata = 1
       seqtype = 0    * 0: nucleotides; 1:codons; 2:AAs
       usedata = 3    * 0: no data; 1:seq like; 2:normal approximation; 3:out.BV (in.BV)
         clock = 2    * 1: global clock; 2: independent rates; 3: correlated rates
       RootAge = '>11<22'  * safe constraint on root age, used if no fossil for root.

         model = 7    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
         alpha = 0.5    * alpha for gamma rates at sites
         ncatG = 4    * No. categories in discrete gamma

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?
