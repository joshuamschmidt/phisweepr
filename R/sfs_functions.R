#' Calculate the 1-d SFS
#'
#' @param genotype_matrix A dgCMat sparse Matrix
#' @param fixed_sites A string
#' @return A list object: A numeric vector of the derived allele counts and the 
#' sfs from  \code{genotype_matrix}
#' @examples
#' get_sfs(genotype_matrix, fixed_sites = FALSE)
get_sfs <- function(genotype_matrix, fixed_derived = FALSE){
  derived_allele_counts <- Matrix::diff(genotype_matrix@p)
  sfs <- table(derived_allele_counts)
  sfs[1] <- 0
  if( fixed_derived == FALSE ){
    sfs[length(sfs)] <- 0
  }
  sfs <- as.vector(sfs/sum(sfs))
  list(derived_allele_counts=derived_allele_counts, sfs=sfs)
}
#' Estimate the mutation rate
#'
#' @param genotype_matrix A list
#' @param fixed_sites A string
#' @return A list object: A numeric vector of the derived allele counts and the 
#' sfs from  \code{genotype_matrix}
#' @examples
#' estimate_mu(genotype_matrix)
estimate_mu <- function(derived_allele_counts,genoObject) {
  theta_from_pi(derived_allele_counts,genoObject$sample.size,genoObject$locusLength)
}

# lapply(0:30, function(x) {sfs <- 1/(0:x); sfs[1] <- 0; sfs})
#----------------
# 2-d SFS functions
# formally called background_jointSFS()
# R version of geting 2-d SFS
#' Estimate the conditional 2d sfs from genotype matrix
#'
#' @param dataObject A list containing a sparse matrix genotypes
#' @param n1Min An int
#' @param n1Max An int
#' @param monomorphic A string
#' @return A data table of n1 k1 k2 n2 counts derived from \code{genotype_matrix}
#' @examples
#' estimate_mu(dataObject, 9, 21, FALSE)

get_two_dimensionalSFS <- function(dataObject,
                                   n1Min,
                                   n1Max,
                                   monomorphic){
  # define output matrix
  # outMatrix <- Matrix::Matrix(0,nrow = dataObject$sample.size+1, ncol =  dataObject$sample.size+1,sparse = F)
  derived_allele_counts <- Matrix::diff(dataObject$genotypes@p)
  # how many sits are due to outgroup (monomorphic but in the genoMatrix)
  n.no.derived <- length(which(derived_allele_counts == 0))
  # therefore the number of monomorphic positions is...
  n.monomorphic <- dataObject$phys_pos[length(dataObject$phys_pos)] - dataObject$phys_pos[1] - length(derived_allele_counts) + n.no.derived
  allDT <- data.table()
  for(i in seq_along(derived_allele_counts)) {
    n1 <- derived_allele_counts[i] 
    if (n1 >= n1Min && n1 <= n1Max) {
      n2 <- dataObject$sample.size - n1
      rowslice <- which(dataObject$genotypes[,i]!=0)
      k1 <- factor(diff(dataObject$genotypes[rowslice,]@p),levels = 0:dataObject$sample.size)
      k2 <- factor(diff(dataObject$genotypes[-rowslice,]@p),levels = 0:dataObject$sample.size)
      # how many positions were monomorphic? add?
      # also realise if simuated with outgroup, can get k1 0 k2 0 positions.
      kDT <- unique(data.table(n1,k1,k2)[!i,][,count:= .N, by=.(n1,k1,k2)])
      kDT <- kDT[!(k1==0 & k2==0)]
      if (monomorphic == TRUE) {
        kDT <- rbindlist(list(kDT,data.table(n1,k1=0,k2=0,count=n.monomorphic)))
      }
      allDT <- rbindlist(list(allDT,kDT))
    }
  }
  allDT <- unique(allDT[,count:= sum(count), by=.(n1,k1,k2)])
  return(allDT)
}
