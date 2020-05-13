ml_phir <- function(core_idx,
                    dataObject,
                    phiNTable,
                    phiSTable,
                    alpha_levels,
                    method="bruteforce",
                    low_alpha,
                    pEscapeLim,
                    log10_alphad.vec,
                    minN1,
                    monomorphic=FALSE) {
  n1 <- dataObject$derived_allele_counts[core_idx]
  n2 <- dataObject$sample.size - n1
  core_pos <- dataObject$phys_pos[core_idx]
  core_window_idxs <- which(abs(dataObject$phys_pos-core_pos) <= distance_pESingle_C(10^low_alpha,1,pEscapeLim) )
  # idx of core within core_window_idxs
  core_idx_in_window <- which(core_window_idxs==core_idx)
  #core_n1 <- test$n1counts[core_snps_idxs[i]]
  rowslice <- which(dataObject$genotypes[,core_idx]!=0)
  k1 <- diff(dataObject$genotypes[rowslice,core_window_idxs]@p)[-core_idx_in_window]
  k2 <- diff(dataObject$genotypes[-rowslice,core_window_idxs]@p)[-core_idx_in_window]
  pos.vec <- dataObject$phys_pos[core_window_idxs][-core_idx_in_window]
  dist.vec <- abs(pos.vec-core_pos)
  if (method=='bruteforce'){
    mle <- maxLRatio_phi_BruteForce(n1,
                                    n2,
                                    k1,
                                    k2,
                                    minN1,
                                    dist.vec,
                                    phiNTable,
                                    phiSTable,
                                    alpha_levels,
                                    pEscapeLim,
                                    log10_alphad.vec,
                                    monomorphic=FALSE,
                                    mc.cores=1,
                                    returnRaw = FALSE)
  }
  mle[,core_pos:= eval(core_pos)]
  return(mle)
}


maxLRatio_phi_BruteForce = function(n1,
                                    n2,
                                    k1,
                                    k2,
                                    minN1,
                                    dist.vec,
                                    phiNTable,
                                    phiSTable,
                                    alpha_levels,
                                    pEscapeLim,
                                    log10_alphad.vec,
                                    monomorphic=FALSE,
                                    mc.cores = 1,
                                    returnRaw = FALSE) {
  optimFun = function(x) -LratioCalcPhi(n1,n2,k1,k2,minN1,phiNTable,phiSTable,pEscapeLim,dist.vec,log10_alphad.vec,monomorphic=FALSE,parVec= c(10^x[1], x[2]))
  optimSol = gridSearch(fun = optimFun, levels = list(seq(alpha_levels[1],alpha_levels[2],length.out = alpha_levels[3]), 1), method = "multicore", mc.control = list(mc.cores = mc.cores),printDetail=FALSE)
  #optimSol$values = log(-optimSol$values)
  optimSol$values = (-1*optimSol$values)
  if (returnRaw == TRUE) {
    return(list(plain = data.table(n1=n1,n2=n2, log10_alpha = optimSol$minlevels[1], beta = optimSol$minlevels[2], LLratio = (-1*optimSol$minfun)), raw = optimSol))
  } else {
    return(data.table(n1=n1,n2=n2, log10_alpha = optimSol$minlevels[1], beta = optimSol$minlevels[2], LLratio = (-1*optimSol$minfun)))
  }
}



LratioCalcPhi = function(n1,
                         n2,
                         k1,
                         k2,
                         minN1,
                         phiNTable,
                         phiSTable,
                         pEscapeLim,
                         dist.vec,
                         log10_alphad.vec,
                         monomorphic=FALSE,
                         parVec = c(alpha, beta)) {
  test_alpha <- round(parVec[1],10)
  test_beta <- parVec[2]
  alpha_core_sub_idxs <- subsetKcountWindow_probEscape_C(dist.vec,test_alpha,test_beta,pEscapeLim)
  # if none of the sites are close enough for this cutoff.
  if (alpha_core_sub_idxs[2] < alpha_core_sub_idxs[1] ){
    return(1) # check this is the right value to return! means 
  }
  if (alpha_core_sub_idxs[2]-alpha_core_sub_idxs[1] < 6) {
    return(1) # check this is the right value to return!
  }
  sub_partition_k1_counts <- k1[alpha_core_sub_idxs[1]:alpha_core_sub_idxs[2]]
  sub_partition_k2_counts <- k2[alpha_core_sub_idxs[1]:alpha_core_sub_idxs[2]]
  sub_dist.vec <- dist.vec[alpha_core_sub_idxs[1]:alpha_core_sub_idxs[2]]
  sub_alphaD.vec = round(log10(sub_dist.vec*test_alpha), 2)
  sub_alphaD.vec[which(sub_alphaD.vec < log10_alphad.vec[2])] <- log10_alphad.vec[1]
  sub_alphaD.vec[which(sub_alphaD.vec > log10_alphad.vec[length(log10_alphad.vec)])] <- log10_alphad.vec[length(log10_alphad.vec)]
  sub_alphaD_idx <- match(sprintf("%0.2f", sub_alphaD.vec), sprintf("%0.2f", log10_alphad.vec))-1
  LLvec_sel <- (log(getProbVectorPhiS(phitable = phiSTable, n1 = n1, minN1 = minN1,k1 = sub_partition_k1_counts,k2 = sub_partition_k2_counts,sub_alphaD_idx = sub_alphaD_idx)))
  LLvec_neu <- (log(getProbVectorPhiN(phitable = phiNTable, n1 = n1, minN1 = minN1,k1 = sub_partition_k1_counts,k2 = sub_partition_k2_counts)))
  if(monomorphic==TRUE){
    centre_idx = which(sub_dist.vec==min(sub_dist.vec))
    mono_lhs = sub_dist.vec[1]:1
    mono_lhs = mono_lhs[!mono_lhs %in% sub_dist.vec[1:centre_idx]]
    mono_rhs = 1:sub_dist.vec[length(sub_dist.vec)]
    mono_rhs = mono_rhs[!mono_rhs %in% sub_dist.vec[(centre_idx+1):length(sub_dist.vec)]]
    mono_dist.vec=c(mono_lhs,mono_rhs)
    mono_alphaD.vec = round(log10(mono_dist.vec*test_alpha), 2)
    mono_alphaD.vec[which(mono_alphaD.vec < log10_alphad.vec[2])] <- log10_alphad.vec[1]
    mono_alphaD.vec[which(mono_alphaD.vec > log10_alphad.vec[length(log10_alphad.vec)])] <- log10_alphad.vec[length(log10_alphad.vec)]
    mono_alphaD_idx <- match(sprintf("%0.2f", mono_alphaD.vec), sprintf("%0.2f", log10_alphad.vec))-1
    LLvec_sel <- LLvec_sel+ sum(log(getProbVectorPhiS_mono(phitable = phiSTable, n1 = n1, minN1 = minN1,sub_alphaD_idx = mono_alphaD_idx)))
    LLvec_neu <- LLvec_neu + sum(getProbVectorPhiN_mono(phitable = obsPhiNTable,n1 = n1,minN1 = minN1,nSites = length(mono_alphaD_idx)))
  }
  # exp(s)
  # CL_sel = exp(LLvec_sel)
  # CL_neu = exp(LLvec_neu)
  # return(2*(log(CL_sel)-log(CL_neu)))
  #return(exp(sum(LLvec_sel) - sum(LLvec_neu)))
  return(sum(LLvec_sel) - sum(LLvec_neu))
}