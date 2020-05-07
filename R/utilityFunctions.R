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

