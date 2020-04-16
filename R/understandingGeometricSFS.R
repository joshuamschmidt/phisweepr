# testing hte calcuation of p[HOMO]
# is actually the summed Pr that a site is polymorophic
# used as 
# if( ktt==0 || ktt==nSam ){
#   LnL1 += log(1.0-theta*Table[n1t][z][0][0]);
#   LnN += log(1.0-theta*pHomo[n1t]); // 1-example.data$mutation_rate*pHomo
# }

nSam <- 30
n1 <- 15
n2 <- nSam - n1
pHomo <- 0
Cb <- binomial.coeffs
for(i in 0:n1) {
  for(j in 0:n2) {
    if( (i+j)!=0 && (i+j)!=nSam ) {
      pHomo <- pHomo + Cb[n1+1,i+1]*Cb[n2+1,j+1]/Cb[nSam+1,i+j+1]/(i+j);
    }
  }
}

pHomon15 <- 0
i=1
for(j in 0:n2) {
  if( (i+j)!=0 && (i+j)!=nSam ) {
    pHomon15 <- pHomon15 + Cb[n1+1,i+1]*Cb[n2+1,j+1]/Cb[nSam+1,i+j+1]/(i+j);
  }
}

pHomon152 <- 0
i=2
for(j in 0:n2) {
  if( (i+j)!=0 && (i+j)!=nSam ) {
    pHomon152 <- pHomon152 + Cb[n1+1,i+1]*Cb[n2+1,j+1]/Cb[nSam+1,i+j+1];
  }
}

# this would be the same as subP[n1] if we also included fixed derived.....
      