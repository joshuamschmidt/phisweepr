# ML config functions
maxd_byAlpha <- function(alphaParVec = c(-6,-2,40), pEscapeLim) {
  log10_alpha_seq = data.table(alpha=10^(seq(alphaParVec[1], alphaParVec[2], length.out = alphaParVec[3])))
  log10_alpha_seq[,max_d:=round(distance_pESingle_C(alpha,1,pEscape = pEscapeLim)),by=alpha]
  return(log10_alpha_seq)
}
