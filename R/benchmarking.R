# sourceCpp("src/armadillo_colSum.cpp")
# 
# library(microbenchmark)
# library(Matrix)
# 
# sp_geno <- example.data$geno
# full_geno <- as.matrix(Matrix(sp_geno,sparse = F))
# microbenchmark(
#   Arma_colSums(full_geno),
#   Sugar_colSums(full_geno),
#   Cpp_colSums(full_geno),
#   Matrix::colSums(full_geno),
#   diff(sp_geno@p),
#   colSums(sp_geno),
#   Matrix::colSums(sp_geno),
#   times=100
# )
# 
# subSparse <- sp_geno[1:10,1:10]
