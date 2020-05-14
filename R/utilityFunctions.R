# Create log scale grid
# from Chris
seq.log10 <- function(from = 1, to = 1, length.out = 50, add.zero = FALSE, shifting = 0, ...){
  res = 10^(seq(from = log10(from + shifting), to = log10(to + shifting), length=length.out - add.zero)) - shifting
  if (add.zero) {
    if (from > to) {
      res = c(res,0)
    } else {
      res = c(0,res)
    }
  }
  res
}

# matrix to longDT
matrixToLongDT <- function(matrixslice,n1,n2){
  longDT <- data.table()
  for(i in 1:(n1+1)){
    for(j in 1:(n2+1)){
      longDT <- rbindlist(list(longDT,data.table(n1=n1,k1=i-1,n2=n2,k2=j-1,pr=matrixslice[j,i])))
    }
  }
  return(longDT)
}


scaled_parmas_sim <- function(mu,r,s,locuslength,Ne){
  s_theta = 4 * Ne * mu * locuslength
  s_rho = 4 * Ne * r * locuslength
  s_sel = 2*Ne*s
  return(c(s_theta,s_rho,s_sel))
  }