#' Calculate the derived allele count per site
#' @param genotype_matrix A dgCMat sparse Matrix
#' @return A numeric vector of the derived allele counts  \code{genotype_matrix}
#' @examples
#' count_derived_alleles(genotype_matrix)
count_derived_alleles <- function(genotype_matrix){
  da_counts <- diff(genotype_matrix@p)
  da_counts[1] <- 0
  return(da_counts)
}


#' Calculate the 1-d SFS
#'
#' @param genotype_matrix A dgCMat
#' @param fixed_sites A string
#' @return A numeric vector of the empirical sfs from  \code{genotype_matrix}
#' @examples
#' get.sfs(genotype_matrix, fixed_sites = FALSE)
get.sfs <- function(genotype_matrix, fixed_derived = FALSE){
  sfs <- table(diff(genotype_matrix@p))
  sfs[1] <- 0
  if( fixed_derived == FALSE){
    sfs[length(sfs)] <- 0
  }
  sfs <- sfs/sum(sfs)
  return(sfs)
}