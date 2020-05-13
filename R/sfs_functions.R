#' Calculate the derived allele count.
#'
#' @param genotype_matrix A dgCMat sparse Matrix
#' @return A vector: A numeric vector of the derived allele counts from \code{genotype_matrix}
#' @examples
#' get_sfs(genotype_matrix)
get_DACount <- function(dataObject, monomorphic=FALSE, fixed_derived = FALSE){
  # max n, use factor to get all counts...
  maxN <- dim(dataObject$genotypes)[1]
  derived_allele_counts <- factor(dataObject$derived_allele_counts,levels = 0:maxN)
  sfs <- table(derived_allele_counts)
  if ( monomorphic == FALSE ){
    sfs[1] <- 0
  }
  if( fixed_derived == FALSE ){
    sfs[length(sfs)] <- 0
  }
  sfs <- as.vector(sfs/sum(sfs))
  return(sfs)
}


#' Calculate the 1-d SFS and derived allele count.
#'
#' @param genotype_matrix A dgCMat sparse Matrix
#' @param fixed_derived A string
#' @param monomorphic A string
#' @return A list object: A numeric vector of the derived allele counts and the 
#' sfs from  \code{genotype_matrix}
#' @examples
#' get_sfs(genotype_matrix, fixed_sites = FALSE)
get_sfsDcountList <- function(dataObject, monomorphic=FALSE, fixedDerived = FALSE){
  # max n, use factor to get all counts...
  maxN <- dim(dataObject$genotypes)[1]
  #derived_allele_counts <- Matrix::diff(genotype_matrix@p)
  derived_allele_counts <- factor(Matrix::diff(dataObject$genotypes@p),levels = 0:maxN)
  sfs <- table(derived_allele_counts)
  if ( monomorphic == FALSE ){
    sfs[1] <- 0
  }
  if( fixedDerived == FALSE ){
    sfs[length(sfs)] <- 0
  }
  if ( monomorphic == TRUE) {
    n.no.derived <- length(which(derived_allele_counts == 0))
    # therefore the number of monomorphic positions is...
    n.monomorphic <- dataObject$phys_pos[length(dataObject$phys_pos)] - dataObject$phys_pos[1] - length(derived_allele_counts) + n.no.derived
    sfs[1] <- n.monomorphic
  }
  sfs <- as.vector(sfs/sum(sfs))
  list(derived_allele_counts=factorToInt(derived_allele_counts), sfs=sfs)
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

get_two_dimensionalSFSDT <- function(dataObject,
                                   n1Min,
                                   n1Max,
                                   monomorphic=FALSE,
                                   fixedDerived=FALSE){
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
      if (n1==1) {
        k1 <- factor(dataObject$genotypes[rowslice,],levels = 0:n1Max)[-i]
      } else{ k1 <- factor(diff(dataObject$genotypes[rowslice,]@p),levels = 0:n1Max)[-i] }
      if (n2==1) {
        k2 <- factor(dataObject$genotypes[-rowslice,],levels = 0:n1Max)[-i]
      } else{ k2 <- factor(diff(dataObject$genotypes[-rowslice,]@p),levels = 0:n1Max)[-i] }
      # how many positions were monomorphic? add?
      # also realise if simuated with outgroup, can get k1 0 k2 0 positions.
      kDT <- unique(data.table(n1,k1,k2)[!i,][,count:= .N, by=.(n1,k1,k2)])
      nFixedDerived <- ifelse(dim(kDT[k1==n1 && k2==dataObject$sample.size-n1])[[1]][1]==0,0L,dim(kDT[k1==n1 && k2==dataObject$sample.size-n1])[[1]][1]==0)
      if (monomorphic == TRUE && fixedDerived== TRUE) {
        kDT <- rbindlist(list(kDT,data.table(n1,k1=0,k2=0,count=n.monomorphic - nFixedDerived)))
      }
      if (monomorphic == TRUE && fixedDerived== FALSE) {
        kDT <- kDT[!(k1==n1 & k2==dataObject$sample.size-n1)]
        kDT <- rbindlist(list(kDT,data.table(n1,k1=0,k2=0,count=n.monomorphic)))
      }
      if (monomorphic == FALSE && fixedDerived== FALSE) {
        kDT <- kDT[!(k1==n1 & k2==dataObject$sample.size-n1)]
        kDT <- kDT[!(k1==0 & k2==0)]
      }
      kDT[,`:=`(k1=factorToInt(k1),k2=factorToInt(k2))]
      if(kDT[,max(as.numeric((k1)))] > n1) {
        print(i)
      }
      allDT <- rbindlist(list(allDT,kDT))
    }
  }
  allDT <- unique(allDT[,count:= sum(count), by=.(n1,k1,k2)])
  allDT <- allDT[,.(k1,k2,pr=count/sum(count)),by=n1][order(n1,k1,k2)]
  return(allDT)
}

get_two_dimensionalSFSlist <- function(dataObject,
                                   n1Min,
                                   n1Max,
                                   monomorphic=FALSE,
                                   fixedDerived=FALSE){
  # this may be unconventional, as I store these as k2,k1 matrices.
  # define output matrix
  # matirx only needs ot cover range 0..n1Max = n1Max+1
  outMatrix <- Matrix::Matrix(0,nrow = n1Max+1, ncol =  n1Max+1,sparse = T)
  # we need n1Max -n1Min + 1
  outList <- rep(list(outMatrix),n1Max - n1Min + 1)
  for(i in seq_along(dataObject$derived_allele_counts)) {
    n1 <- dataObject$derived_allele_counts[i]
    if (n1 >= n1Min && n1 <= n1Max) {
      n <- n1 - n1Min + 1 # relative n, to find which list slice stores the n1 matrix.
      n2 <- dataObject$sample.size - n1
      rowslice <- which(dataObject$genotypes[,i]!=0)
      # this breaks down for n=1!!!!!
      if (n1==1) {
        k1 <- factor(dataObject$genotypes[rowslice,],levels = 0:n1Max)[-i]
      } else{ k1 <- factor(diff(dataObject$genotypes[rowslice,]@p),levels = 0:n1Max)[-i] }
      if (n2==1) {
        k2 <- factor(dataObject$genotypes[-rowslice,],levels = 0:n1Max)[-i]
      } else{ k2 <- factor(diff(dataObject$genotypes[-rowslice,]@p),levels = 0:n1Max)[-i] }
      subSFS <- unclass(table(k2,k1))
      if(monomorphic == TRUE && fixedDerived== TRUE) {
        n.no.derived <- subSFS[1,1]
        n.derived <- sum(subSFS) - subSFS[1,1]
        # therefore the number of monomorphic positions is...
        n.monomorphic <- dataObject$phys_pos[length(dataObject$phys_pos)] - dataObject$phys_pos[1] - n.derived + n.no.derived
        subSFS[1,1] <- n.monomorphic
      }
      if(monomorphic == TRUE && fixedDerived== FALSE) {
        n.no.derived <- subSFS[1,1]
        n.f.derived <- subSFS[n2+1,n1+1]
        subSFS[n2+1,n1+1] <- 0
        n.derived <- sum(subSFS) - subSFS[1,1]
        n.monomorphic <- dataObject$phys_pos[length(dataObject$phys_pos)] - dataObject$phys_pos[1] - n.derived + n.no.derived - n.f.derived
        subSFS[1,1] <- n.monomorphic
      }
      if(monomorphic == FALSE && fixedDerived== FALSE) {
        subSFS[n2+1,n1+1] <- 0
        subSFS[1,1] <- 0
      }
      outList[[n]] <- outList[[n]] + subSFS
    }
    
  }
  outList <- lapply(outList, function(x) x/sum(x))
  if(fixedDerived== FALSE){
    outList[[dataObject$sample.size]] <- Matrix::Matrix(0,nrow = n1Max+1, ncol =  n1Max+1,sparse = T)
    outList[[dataObject$sample.size]][1,] <- dataObject$sfs
  }
  return(outList)
}


getSfs_perN_from_two_dimensionalSFSlist <- function(two_dimensionalSFSlist,
                                                    monomorphic=FALSE,
                                                    fixedDerived=FALSE){
  list_neutral_sfs <- mapply(calc_MatrixCell_divededbyColSum, two_dimensionalSFSlist, seq_along(two_dimensionalSFSlist),MoreArgs=list(monomorphic,fixedDerived),SIMPLIFY = FALSE )
  names(list_neutral_sfs) = 1:length(list_neutral_sfs)
  list_neutral_sfs <- append(list_neutral_sfs, list('0'=0), 0)
  sfs_Matrix <- matrix(0.0,nrow=length(list_neutral_sfs),ncol=length(list_neutral_sfs))
  for(i in 1:length(list_neutral_sfs)) {
    # so R 1 is cpp 0/SFS n0
    n_sfs <- list_neutral_sfs[[i]]
    if(length(n_sfs) >0) {
      sfs_Matrix[i, 1:length(n_sfs)] <- n_sfs
    }
  }
  return(sfs_Matrix) 
}

calc_MatrixCell_divededbyColSum <- function(two_dimensionalSFS, n, monomorphic= FALSE, fixedDerived= FALSE){
  sfs <-  Matrix::colSums(two_dimensionalSFS)
   if(monomorphic == FALSE) {
    sfs[1] <- 0
   }
  if(fixedDerived == FALSE) {
    sfs[n+1] <- 0
  }
  if(sum(sfs)> 0){
    sfs <- sfs/sum(sfs)
  }
  sfs <- sfs[1:(n+1)]
  return(sfs)
  }


factorToInt <- function(intFactorVector){
  factorAsIntVec <- as.integer(levels(intFactorVector))
  converted <- factorAsIntVec[intFactorVector]
  return(converted)  
}

sfsMatrixToLongDT <- function(sfsMatrix, n1, nSam) {
  # construct output DT
  n2 <- nSam - n1
  outDT <- data.table(n1=eval(n1),expand.grid(k1=0:n1,k2=0:n2))
  outDT[,frequency:= 0]
  if(is.data.table(outDT)==FALSE){
    stop("not a data.table!")
  } else{
    print("a data.table!")
  }
  # fill DT with values from Matrix
  for(i in 1:(n1+1)){
    for(j in 1:(n2+1)){
      outDT[k1 == (i-1) & k2 == (j-1)]$frequency <- sfsMatrix[j, i]
    }
  }
  return(outDT)
} 

filterDataObject <- function(dataObject, monomorphic= FALSE, fixedDerived= FALSE) {
  if(monomorphic == FALSE && fixedDerived == FALSE) {
    unfiltered_idxs <- which(dataObject$derived_allele_counts!=0 & dataObject$derived_allele_counts!=dataObject$sample.size)
  }
  if(monomorphic == FALSE && fixedDerived == TRUE) {
    unfiltered_idxs <- which(dataObject$derived_allele_counts!=0)
  }
  if(fixedDerived == FALSE && monomorphic == TRUE) {
    unfiltered_idxs <- which(dataObject$derived_allele_counts!=dataObject$sample.size)
  }
  if(fixedDerived == TRUE && monomorphic == TRUE) {
    return(dataObject)
  }
  dataObject$phys_pos <- dataObject$phys_pos[unfiltered_idxs]
  dataObject$gen_pos <- dataObject$gen_pos[unfiltered_idxs]
  dataObject$genotypes <- dataObject$genotypes[,unfiltered_idxs]
  dataObject$derived_allele_counts <- dataObject$derived_allele_counts[unfiltered_idxs]
  return(dataObject)
}





getSfs_perN_from_downprojectFullSFS <- function(sfs,scaled=TRUE,monomorphic= FALSE, fixedDerived= FALSE){
  n = length(sfs) - 1
  list_neutral_sfs = lapply(1:n, function(x) {sfs_raw = p_jH_vec_C(pvec = sfs, j = 0:x, H = x, n = n)})
  #names(list_neutral_sfs) = 1:n
  if(monomorphic == FALSE && fixedDerived == FALSE) {
    list_neutral_sfs <- mapply(function(x, i)  {x[1]<-0;x[i+1]<-0; return(x);}, list_neutral_sfs, seq_along(list_neutral_sfs),SIMPLIFY = FALSE )
    }
  if(monomorphic == FALSE && fixedDerived == TRUE) {
    list_neutral_sfs <- mapply(function(x, i)  {x[1]<-0; return(x);}, list_neutral_sfs, seq_along(list_neutral_sfs),SIMPLIFY = FALSE )
  }
  if(fixedDerived == FALSE && monomorphic == TRUE) {
    list_neutral_sfs <- mapply(function(x, i)  {x[i+1]<-0; return(x);}, list_neutral_sfs, seq_along(list_neutral_sfs),SIMPLIFY = FALSE )
  }
  list_neutral_sfs <- append(list_neutral_sfs, list('0'=0), 0)
  sfs_Matrix <- matrix(0.0,nrow=length(sfs),ncol=length(sfs))
  for(i in 1:length(list_neutral_sfs)) {
    # so R 1 is cpp 0/SFS n0
    n_sfs <- list_neutral_sfs[[i]]
    if(scaled==TRUE & i > 2){
      n_sfs <- n_sfs/sum(n_sfs)
    }
    sfs_Matrix[i, 1:length(n_sfs)] <- n_sfs
  }
  return(sfs_Matrix)
}


getSfs_perN_from_downsampleFullSFS <- function(dataObject,monomorphic= FALSE, fixedDerived= FALSE){
  n = dim(dataObject$genotypes)[1]
  nvec=1:n
  list_neutral_sfs = mapply(sfsResampleGenotypeMatrix, i=nvec, nSam=n,MoreArgs=list(genotypes=dataObject$genotypes,phys_pos=dataObject$phys_pos,monomorphic,fixedDerived),SIMPLIFY = FALSE )
  list_neutral_sfs <- append(list_neutral_sfs, list('0'=0), 0)
  sfs_Matrix <- matrix(0.0,nrow=n+1,ncol=n+1)
  for(i in 1:length(list_neutral_sfs)) {
    # so R 1 is cpp 0/SFS n0
    n_sfs <- list_neutral_sfs[[i]]
    # if(scaled==TRUE & i > 2){
    #   n_sfs <- n_sfs/sum(n_sfs)
    # }
    sfs_Matrix[i, 1:length(n_sfs)] <- n_sfs
  }
  return(sfs_Matrix)
}


sfsResampleGenotypeMatrix <- function(genotypes,phys_pos,i, nSam,monomorphic=FALSE,fixedDerived=FALSE){
  rowslice <- sort(sample(1:nSam,size = i,replace = FALSE))
  if(i==1){
    sfs=table(genotypes[rowslice,])
  } else{sfs=table(factor(diff(genotypes[rowslice,]@p),levels = 0:i))}
  # check that the last value, i+1 exists
  if(monomorphic==FALSE){
    sfs[1] <- 0
  }
  if(fixedDerived==FALSE){
    sfs[i+1] <- 0
  }
  if(monomorphic ==TRUE) {
    n.no.derived <- sfs[1]
    n.derived <- sum(sfs[-1])
    # therefore the number of monomorphic positions is...
    n.monomorphic <-phys_pos[length(phys_pos)] - phys_pos[1] - n.derived + n.no.derived
    sfs[1] <- n.monomorphic
  }
  if(i>1){
    sfs <- sfs/sum(sfs)
  }
  return(sfs)
}

make_log10alphad <- function(minalpha,maxalpha,length.out,add.zero = TRUE){
  alphad_seq = seq.log10(minalpha, maxalpha, length.out = length.out, add.zero = TRUE)
  log10_alphad=round(log10(alphad_seq),2)
  return(log10_alphad)
}

alphaDres*100+2

getPhiS <- function(dataObject, nkMatrix, n1range, minalpha,maxalpha,length.out,add.zero = TRUE, alphaDres=alphaDres,log10_alphad=,beta=1){
  
  phiSTable <- makePhiSTable(nSam = dataObject$sample.size,testN1s = n1range, ptable = nkMatrixds, alphad = 10^log10_alphad,beta)
  return(phiSTable)
}
