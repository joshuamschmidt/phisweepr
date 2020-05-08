# ML functions

LratioCalcPhi = function(k1counts,
                         k2counts,
                         dist.vec,
                         n1,
                         n2,
                         phiNtable,
                         phiDtable,
                         alphad_comb_list,
                         mono_sites=FALSE,
                         parVec = c(alpha, beta,pEscapeLim)) {
  test_alpha <- round(parVec[1],10)
  test_beta <- parVec[2]
  alpha_core_sub_idxs <- subsetKcountWindow_probEscape_C(partition_k_counts$d,
                                                         test_alpha,
                                                         test_beta,
                                                         parVec[3])
  # Only 5 SNPs, plus core SNP, is not enough
  #if (length(alpha_core_sub_idxs) < 6) {
  #    return(1)
  #}
  if (alpha_core_sub_idxs[2]-alpha_core_sub_idxs[1] < 6) {
    return(1)
  }
  sub_partition_k_counts <- partition_k_counts[alpha_core_sub_idxs[1]:alpha_core_sub_idxs[2],]
  sub_partition_k_counts <- addMonomorphicPhiS(sub_partition_k_counts,
                                               mono_sites)
  # fold the distances
  maxd <- alphad_comb_list[[1]][alpha==eval(test_alpha),max_d]
  sub_partition_k_counts[d > maxd,d:=maxd]
  LLvec_neu <- log(phiN_hash[[sub_partition_k_counts[,paste(eval(n1),k1,
                                                            eval(n2),k2)]]])
  sub_partition_k_counts<- sub_partition_k_counts[d!=0,log10_alphad:=round(alphad_comb_list[[2]][[paste(test_beta,test_alpha,d)]],2) ]
  system.time(LLvec_sel <- log(phiS_hash[[sub_partition_k_counts[d!=0,paste(n1,k1,
                                                                            n2,k2,
                                                                            log10_alphad,
                                                                            test_beta)]]]))
  
  system.time(LLvec_sel_f <- log(phiS_hash[[sub_partition_k_counts[d!=0,paste(n1,k1,
                                                                              n2,k2,
                                                                              round(get_log10alphad_C(test_alpha,d),2),
                                                                              test_beta)]]]))
  return(exp(sum(LLvec_sel) - sum(LLvec_neu)))
}