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
get_jointSFS <- function(dataObject, testn1) {
  n1_idxs <- which(dataObject$derived_allele_counts == testn1)
  if(length(n1_idxs)==0) {
    return(NA)
  }
  n2 <- dataObject$sample.size - testn1
  # idx positions for sites at n1 frequency
  partition_k_counts <- data.table()
  for (i in 1:length(n1_idxs)) {
    idx <- n1_idxs[i]
    if (i==1) {
      partition_k_counts <- data.table(groupInfo(genoMatrix,idx)[-idx,])[,.N,by=.(V1,V2)][order(V1,V2)]
    }
    else {
      partition_k_counts <- rbindlist(list(partition_k_counts, data.table(groupInfo(genoMatrix,idx)[-idx,])[,.N,by=.(V1,V2)][order(V1,V2)]))
    }
  }
  names(partition_k_counts) <- c("k1","k2","N")
  partition_k_counts <- partition_k_counts[,.(N=sum(N)),by=.(k1,k2)][order(k1,k2)]
  #missing entries in the joint sfs?
  setkey(partition_k_counts,k1,k2)
  all_possibleSFS <- data.table(expand.grid(x= 0:testn1, y=0:n2))[,.(k1=x,k2=y)][order(k1,k2)]
  setkey(all_possibleSFS,k1,k2)
  partition_k_counts <- partition_k_counts[all_possibleSFS]
  partition_k_counts[N %in% NA,N:=0]
  partition_k_counts[,`:=`(n1=eval(testn1),n2=eval(n2))]
  # scale
  partition_k_counts[,pr:=N/sum(N)]
  # standardising.
  # we have looped over x = length(testn1) positions.
  # thus the total number of variants in the 2-d SFS is # polymorphic
  # * x. so actually unique poly count is partition_k_counts[,sum(N)]/length(n1_idxs).
  # so we multiply proportions to partition_k_counts[,sum(N)]/length(n1_idxs)/locusLength * pr
  partition_k_counts[,scaledpr:=pr*(partition_k_counts[,sum(N)]/length(n1_idxs)/locusLength)]
  # scaledpr for k1=0 and k2=0, will  be same lookupo for k1=n1 and k2=n2
  ks0 <- 1-partition_k_counts[!(k1==0 & k2==0),sum(scaledpr)]
  partition_k_counts[(k1==0 & k2==0),scaledpr:=ks0]
  return(partition_k_counts)
}