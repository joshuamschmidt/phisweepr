# phisweepr


This is an R/Rcpp implementation of Vy and Kim 2015.

In R:
```R
options(scipen=999)
require(pacman)
pacman::p_load(
devtools,
usethis,
roxygen2,
data.table,
ggplot2,
viridis,
cowplot,
Matrix,
NMOF,
parallel,
Rcpp,
RcppArmadillo)
source('R/simdata_to_Robjects.R')
source('R/sfs_functions.R')
source('R/utilityFunctions.R')
source('R/MLconfig_functions.R')
source('R/MLfunctions.R')
Rcpp::sourceCpp('src/rcpp_sfs_functions.cpp')
Rcpp::sourceCpp('src/armaCubeFields.cpp')
Rcpp::sourceCpp('src/MLconfig_cpp_functions.cpp')
Rcpp::sourceCpp('src/ML_compute_cpp.cpp')


# some data, more files in 'inst/extdata'
dataObject <- get.simulated.data("inst/extdata/msms.selection.100chr.s05.b07_50mb.txt", "msms")

#
monomorphic=FALSE
fixedDerived=FALSE
dataObject <- append(dataObject, get_sfsDcountList(dataObject, monomorphic=monomorphic , fixedDerived = fixedDerived))
if(monomorphic==FALSE | fixedDerived==FALSE) {
  dataObject <- filterDataObject(dataObject,monomorphic=monomorphic , fixedDerived = fixedDerived)
}
dataObject[["mutation_rate"]] <- theta_from_pi(dataObject$derived_allele_counts,nSam = dataObject$sample.size, locusLength = dataObject$locusLength)
n1Min <- round(dataObject$sample.size*0.20)
n1Max <- round(dataObject$sample.size*0.80)
n1range <- n1Min:n1Max


# these is how Vy and Kim derived there std neutral exp sfs. more for my personal benefit. though i do use binomial.coeffs in the expSFS calc.
binomial.coeffs <- CbTable(dataObject$sample.size)
subP <- subPop_freqSpec(binomial.coeffs,dataObject$sample.size)
pHomo <- calc_pHomo(binomial.coeffs,dataObject$sample.size,n1Max,n1Min)
## can get the same thing as.....
subP_2 <- Matrix::Matrix(do.call(rbind,lapply(0:dataObject$sample.size,function(x) {sfs <- numeric(dataObject$sample.size+1); if(x > 1) {sfs[1:x] <-  1/(0:(x-1))}; sfs[1] <- 0; sfs })),sparse = TRUE)
# these are very close in values, down to 2e-7 at WORST.
unlist(lapply(1:(dataObject$sample.size+1), function(x) all.equal(subP[x,],subP_2[x,],tolerance = 2e-7)))


## conditional 2d-SFS
# from observed
obs.cond.2d.sfsL <- get_two_dimensionalSFSlist(dataObject = dataObject,n1Min = 1,n1Max = dataObject$sample.size,monomorphic = monomorphic,fixedDerived = fixedDerived) # this function is SLOW, and I do not like it....
# one needs to determine the 2dSFS, n1,k1,n2,k2. As is, this iterates over all polymorphic sites by n1 class, than suims up the 2dSFS over n1 sites. 

# so this does not scale well with increased chromosome length. I am working on using the downsampled or the hypergeometric projected 1d-SFS, but haven't quite solved it yet.

# calculate phiN
# could try rewritting as mapply or mcmapply process.
# arma::cube is returned to R as a 3D-array. subset[n1,k1,n1matrix]
obsPhiNTable <- listMatricesToArmaCube(sfslist = obs.cond.2d.sfsL,n1Range = n1range,nSam = dataObject$sample.size)
# assumed equilibirum sfs
expPhiNTable <- standardNeutralJointSFSArmaCube(CbTable = binomial.coeffs,pHomo = pHomo, n1Min = n1Min, n1Max = n1Max, nSam = dataObject$sample.size,theta = dataObject$mutation_rate, monomorphic = monomorphic,fixedDerived = fixedDerived)

# calculate phiS
# a matrix relating derived allele n, and k i.e rows are n=0,...n=samplesize, cols are k=0...k=n1
# three ways to estimate this.
# as above, condition sites on n1 count, and then see SFS of given N. summed across all sites of n
nkMatrix <- getSfs_perN_from_two_dimensionalSFSlist(two_dimensionalSFSlist = obs.cond.2d.sfsL,monomorphic = monomorphic,fixedDerived = fixedDerived)
# downproject from the full sample size SFS
nkMatrixdp <- getSfs_perN_from_downprojectFullSFS(sfs = dataObject$sfs,scaled = TRUE,monomorphic = monomorphic,fixedDerived = fixedDerived)
# down sample from the full sample size genotype array
nkMatrixds <- getSfs_perN_from_downsampleFullSFS(dataObject = dataObject,monomorphic = monomorphic,fixedDerived = fixedDerived)

alphaDres <- 4 # a control switch to set n. of alphaD gridpoints, 1,2 or 4... 102,202,402 length of log10_alphad
# note still pass param beta, even though vy and kim define model with beta=1. this is to allow easy upgrade of this if we want.
log10_alphad <- make_log10alphad(minalpha=0.001,maxalpha=10,length.out = alphaDres*100+2)
phiSTable <- makePhiSTable(nSam = dataObject$sample.size,testN1s = n1range,ptable = nkMatrixds,alphad = 10^log10_alphad,beta = 1)


# define params for ML calculation
# first is resolution of alpha maximisiation
alpha_levels <- c(-6,-2,40)
pEscapeLim = 0.95 # % pescape of a single lineage


# define sites to sweep over  
core_snps_idxs <- which(dataObject$derived_allele_counts %in% n1range)
low_alpha <- -6


mc_maxLRatio_table <- mcmapply(ml_phir,core_snps_idxs, MoreArgs = list(dataObject = dataObject,phiNTable = expPhiNTable,phiSTable = phiSTable,alpha_levels=alpha_levels,pEscapeLim = pEscapeLim,method="bruteforce",log10_alphad.vec = log10_alphad,low_alpha=low_alpha,monomorphic = monomorphic,minN1 = n1range[1]),SIMPLIFY = FALSE,mc.cores = 4)

maxLRatioDT <- data.table(matrix(unlist(mc_maxLRatio_table), nrow=length(mc_maxLRatio_table), byrow=T))

# still some issues with returning LL. so for now i maximise diffSumLogs: sum(log(sel))-sum(log(neu)) 

names(maxLRatioDT) <- c("n1","n2","alpha","beta","diffSumLogs","position")
maxLRatioDT[,site.class:=ifelse(position==25e6,"sweep","nonsweep")]
size.mod <- c(5,0.75)
alpha.mod <- c(1,0.5)
names(size.mod) <- c("sweep","nonsweep")
names(alpha.mod) <- c("sweep","nonsweep")
ggplot(maxLRatioDT, aes(x=position,y=diffSumLogs,color = site.class)) + geom_point(alpha=alpha.mod[maxLRatioDT$site.class],size=size.mod[maxLRatioDT$site.class]) + theme_bw()
maxLRatioDT[,std:=(diffSumLogs-mean(diffSumLogs))/sd(diffSumLogs),by=n1]
maxLRatioDT[,.(n=.N,m.diffSumLogs=mean(diffSumLogs)),by=n1][order(n1)]
ggplot(maxLRatioDT, aes(x=position,y=std,color = site.class)) + geom_point(alpha=alpha.mod[maxLRatioDT$site.class],size=size.mod[maxLRatioDT$site.class]) + theme_bw()
saveRDS(maxLRatioDT,file='maxLRatioDT.msms.selection.100chr.s05.b07_50mb.rds')

```
